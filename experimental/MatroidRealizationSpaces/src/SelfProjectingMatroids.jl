@attributes mutable struct MatroidRealizationSpaceSelfProjecting{BaseRingType, RingType} <: AbsAffineScheme{BaseRingType, RingType}
  defining_ideal::Union{Ideal,NumFieldOrderIdeal}
  inequations::Vector{RingElem}
  ambient_ring::Ring
  selfprojecting_realization_matrix::Union{MatElem,Nothing}
  char::Union{Int,Nothing}
  q::Union{Int,Nothing}
  ground_ring::Ring
  one_realization::Bool

  # Fields for caching
  underlying_scheme::AbsAffineScheme{BaseRingType, RingType}

  function MatroidRealizationSpaceSelfProjecting(
    I::Union{Ideal,NumFieldOrderIdeal},
    ineqs::Vector{<:RingElem},
    R::Ring,
    mat::Union{MatElem,Nothing},
    char::Union{Int,Nothing},
    q::Union{Int,Nothing},
    ground_ring::Ring)
    BaseRingType = typeof(ground_ring)
    RingType = typeof(R)
    if char !=0
      error("Code is not implemented for characteristic not equal to zero.")
    end
    if q != nothing
      error("Code is not implemented for finite fields.")
    end
    if R isa MPolyRing
      PolyRingType = typeof(R)
      MultSetType = MPolyPowersOfElement{BaseRingType, elem_type(BaseRingType), PolyRingType, elem_type(PolyRingType)}
      RingType = MPolyQuoLocRing{BaseRingType, elem_type(BaseRingType), PolyRingType, elem_type(PolyRingType), MultSetType}
    end
    return new{BaseRingType, RingType}(I, ineqs, R, mat, char, q, ground_ring, false)
  end
end

function underlying_scheme(RS::MatroidRealizationSpaceSelfProjecting{BRT, RT}) where {BRT<:Ring, RT<:MPolyQuoLocRing}
  isdefined(RS, :underlying_scheme) && return RS.underlying_scheme::AffineScheme{BRT, RT}

  P = ambient_ring(RS)::MPolyRing
  I = defining_ideal(RS)::MPolyIdeal
  U = MPolyPowersOfElement(P, P.(inequations(RS)))::MPolyPowersOfElement
  RS.underlying_scheme = spec(P, I, U)
  return RS.underlying_scheme::AffineScheme{BRT, RT}
end

function underlying_scheme(RS::MatroidRealizationSpaceSelfProjecting{BRT, RT}) where {BRT<:Ring, RT<:MPolyQuoRing}
  isdefined(RS, :underlying_scheme) && return RS.underlying_scheme::AffineScheme{BRT, RT}

  RS.underlying_scheme = spec(ambient_ring(RS)) # for some reason this is an MPolyQuoRing already
  return RS.underlying_scheme::AffineScheme{BRT, RT}
end

#Currently there is no analogue for is_realizable and one_realization for MatroidRealizationSpaceSelfProjecting
function Base.show(io::IO, ::MIME"text/plain", RS::MatroidRealizationSpaceSelfProjecting)
  if isone(defining_ideal(RS))
    print(io, "The matroid does not have a self-projecting realization over characteristic zero.")
  else 
    io = Oscar.pretty(io)
    println(io, "The selfprojecting realization space is")
    print(io, Oscar.Indent())
    show(io, MIME("text/plain"), RS.selfprojecting_realization_matrix)
    print(io, "\n", Oscar.Dedent(), "in the ", Oscar.Lowercase(), RS.ambient_ring)
    I = RS.defining_ideal
    if !iszero(I)
      print(io, "\nwithin the vanishing set of the ideal\n", I)
    end
    if length(RS.inequations) > 0
      print(io, "\navoiding the zero loci of the polynomials\n", RS.inequations)
    end
  end
end
@doc raw"""
    defining_ideal(RS::MatroidRealizationSpaceSelfProjecting)

The ideal of the matroid realization space `RS`.
"""
defining_ideal(RS::MatroidRealizationSpaceSelfProjecting) = RS.defining_ideal

@doc raw"""
    inequations(RS::MatroidRealizationSpaceSelfProjecting)

Generators of the localizing semigroup of `RS`.
These are the polynomials that need to be nonzero in any selfprojecting realization.
"""
inequations(RS::MatroidRealizationSpaceSelfProjecting) = RS.inequations

@doc raw"""
    ambient_ring(RS::MatroidRealizationSpaceSelfProjecting)

The polynomial ring containing the ideal `defining_ideal_sp(RS)` and the polynomials in `inequations(RS)`.
"""
ambient_ring(RS::MatroidRealizationSpaceSelfProjecting) = RS.ambient_ring



@doc raw"""
    selfprojecting_realization_matrix(RS::MatroidRealizationSpaceSelfProjecting)

A matrix with entries in ambient_ring_sp(RS) whose columns, when filled in with values satisfying equalities
from `defining_ideal(RS)` and inequations from `inequations(RS)`, form a self-projecting realization for the matroid.
"""
selfprojecting_realization_matrix(RS::MatroidRealizationSpaceSelfProjecting) = RS.selfprojecting_realization_matrix

function satisfies_disjointbasisproperty(mat::Matroid)::Bool
    dmat = dual_matroid(mat)
    boo = false
    for b in bases(dmat)
        boo = false
        for e in bases(mat)
            if issubset(e,b)
                boo = true
                break
            end
        end
        if !boo 
            break
        end
    end
    return boo;
end


@doc raw"""
    is_selfprojecting(mat::Matroid)

Returns a boolean which states whether the given matroid satisfies the property to be self-projecting.
# Examples
```jldoctest
julia> m = fano_matroid()
Matroid of rank 3 on 7 elements

julia> is_selfprojecting(m)
true
```
"""
#A function to test if a matroid is self-projecting
function is_selfprojecting(mat::Matroid)::Bool
  k = rank(mat)
  n = length(matroid_groundset(mat))
  boo_basis = false
  for b in bases(mat)
    boo_basis = false
    for c in bases(mat)
      if intersect(b,c) == []
        boo_basis = true
        break  
      end
    end
    if !boo_basis 
      break
    end
  end
  boo_flats = true
  for f in flats(mat,k-1)
    boo_flats = true
    for g in flats(mat,k-1)
      if length(union(f,g)) == n-1
        boo_flats = false
        break
      end
    end
    if !boo_flats 
      break
    end
  end
  return boo_basis && boo_flats;
end

@doc raw"""
    selfprojecting_realization_ideal(m::Matroid)

  Function to compute the defining_ideal of a selfprojecting realization space

  !!! warning "This function is slow except for small matroids!"
# Examples
```jldoctest
julia> m = matroid_from_nonbases([[1,2,3],[4,5,6]],6)
Matroid of rank 3 on 6 elements

julia> selfprojecting_realization_ideal(m)
Ideal generated by
  0

julia> m = uniform_matroid(3,6)
Matroid of rank 3 on 6 elements

julia> selfprojecting_realization_ideal(m)
Ideal generated by
  x1*x2*x3 - x1*x2*x4 - x1*x3*x4 + x1*x4 + x2*x3*x4 - x2*x3
```
"""
function selfprojecting_realization_ideal(m::Matroid; saturate::Bool = false, check::Bool = true)::Ideal
  if check 
    @req is_selfprojecting(m) "The given matroid is not self projecting" 
  end
  if !is_realizable(m,char = 0)
    return defining_ideal(realization_space(m,char =0))
  end
  RS = realization_space(m,char =0,simplify = true, saturate = true); 
  R = ambient_ring(RS); 
  I = defining_ideal(RS);
  n = length(matroid_groundset(m))
  k = rank(m)
  RR, x,l = polynomial_ring(QQ, :x=> 1:length(gens(R)), :l=>1:n; cached=false)
  F = hom(R, RR,x) #use this to move ideals into RR
  G = hom(RR, R, vcat(gens(R),[0 for i in 1:n])) #use this to move them back again into the ambient ring of RS #Careful! this sets all occurences of l to zero! Use only when l is eliminated !
  L = ideal(RR, prod(l[i] for i in 1:n))
  MRS = realization_matrix(RS)
  M = matrix(RR,nrows(MRS),ncols(MRS),[F(MRS[i,j]) for i in 1:nrows(MRS) for j in 1:ncols(MRS)])
  DL = diagonal_matrix(RR,[l[i] for i in 1:n])
  V = M*DL*transpose(M)
  IS = ideal(RR,[V[i,j] for i in 1:k for j in 1:k]) + ideal(RR,[F(g) for g in gens(I)])
  J = stepwise_saturation(IS,l) 
  E = eliminate(J,l) 
  if saturate 
    E = stepwise_saturation(E,F.(inequations(RS)))
  end
  if is_zero(E)
    return ideal(ambient_ring(RS),[0]);
  elseif is_one(E)
    return ideal(ambient_ring(RS),[1]);
  else
    return ideal(R,G.(gens(E)))
  end
end


#this function sorts the list of polynomials by total degree, this helps later to reduce the list of basis_minors unfortunately no sort_by_degree function is currently available for quotienring elements
function sort_by_degree(polys::Vector{<:MPolyRingElem})::Vector{<:MPolyRingElem}
  deg(p) = total_degree(p) 
  return sort(polys, by=deg)
end
function sort_by_degree(polys::Vector{<:PolyRingElem})::Vector{<:PolyRingElem}
  deg(p) = degree(p) 
  return sort(polys, by=deg)
end




#this function computes the minors of a given matrix for the bases in a given list.
#this list is not reduced, i.e. there are no repetitions but the polynomials are not reduced, i.e. the list might contain x+1 and y-2 and their product
function basis_minors(M::MatElem, Bases::Vector{Vector{Int}})::Vector{<:RingElem}
  R = base_ring(M);
  if R isa MPolyQuoRing
    candidates = [det(M[1:nrows(M),b]) for b in Bases]
  else
    candidates = sort_by_degree([det(M[1:nrows(M),b]) for b in Bases])
  end
  multiplicativeSet= nothing
  ineqs = [R(0)]
  if R isa MPolyQuoRing #in this case one cannot make the multiplicativeSet using the powers_of_element, so currently the simpler choice is used and there might be double entries or products of polynomials in the list - is there a better solution for this?
    for candidate in candidates
      @req !iszero(candidate) "A basis has vanishing minor. Please check that the input vector of bases really consists of bases of the matroid given by the columns of the input matrix."
      if isone(candidate) || isone(-candidate) 
      elseif !(candidate in ineqs[2:end])&& !(-candidate in ineqs[2:end])
        push!(ineqs,R(candidate))
      end
    end
    return ineqs[2:end]
  end
  for candidate in candidates
    @req !iszero(candidate) "A basis has vanishing minor. Please check that the input vector of bases really consists of bases of the matroid given by the columns of the input matrix."
    if isone(candidate) || isone(-candidate) 
    else
      if isnothing(multiplicativeSet) && !(candidate in ineqs[2:end]) && !(-candidate in ineqs[2:end])
        multiplicativeSet = powers_of_element(candidate)
        push!(ineqs,R(candidate))
      elseif !(candidate in ineqs[2:length(ineqs)])&& !(-candidate in ineqs[2:length(ineqs)]) && !(candidate in multiplicativeSet)  && !(-candidate in multiplicativeSet)
        multiplicativeSet = product(multiplicativeSet,powers_of_element(candidate));
        push!(ineqs,R(candidate))
      end
    end
  end
  return ineqs[2:length(ineqs)]
end


@doc raw"""
  selfprojecting_realization_matrix(m::Matroid, Bas::Vector{Int}, F::Ring;Ideal::Union{Ideal,Nothing} = nothing, check::Bool = true)

  Function to compute a template for a selfprojecting realization. Returns the ambient ring and the matrix. The input Bas defines the basis of the matriod that is chosen to be represented by the identity matrix in the realization. 
  The defining_ideal of the selfprojecting realization space of the given matroid can optionally be entered. Otherwise it will be computed. 
  The kwarg check states whether the matroid will be checked for self-projectivity. The default is true.

  !!! warning "This function uses the computation of the selfprojecting_realization_ideal. Therefore, it is slow except for small matroids."
# Examples
```jldoctest
julia> m = matroid_from_nonbases([[1,2,3],[4,5,6]],6)
Matroid of rank 3 on 6 elements

julia> R, M =selfprojecting_realization_matrix(m,[1,2,4])
(Multivariate polynomial ring in 2 variables over QQ, [1 0 1 0 1 1; 0 1 1 0 x1 x1; 0 0 0 1 1 x2])

julia> M
[1   0   1   0    1    1]
[0   1   1   0   x1   x1]
[0   0   0   1    1   x2]

julia> m = uniform_matroid(3,6)
Matroid of rank 3 on 6 elements

julia> R, M =selfprojecting_realization_matrix(m,[1,2,3])
(Quotient of multivariate polynomial ring by ideal (x1*x2*x3 - x1*x2*x4 - x1*x3*x4 + x1*x4 + x2*x3*x4 - x2*x3), [1 0 0 1 1 1; 0 1 0 1 x1 x3; 0 0 1 1 x2 x4])
```
"""
###################
#Bas are the given columns that will be the identity matrix
#This function returns the ambient ring and the realization matrix, same as for the classical function realization_space_matrix
function selfprojecting_realization_matrix(m::Matroid, Bas::Vector{Int}; I::Union{Ideal,Nothing} = nothing, check::Bool = true)
  if check 
    @req is_selfprojecting(m) "The given matroid is not self projecting" 
  end
  if !(is_realizable(m,char=0))
    return (nothing, nothing)
  else
  X = realization_matrix(realization_space(m,B=Bas,ground_ring=QQ)) #this matrix is not yet simplified by the defining ideal of the selfprojecting realization ideal
    if isnothing(I) #this way the ideal needs only be computed once!
      I = selfprojecting_realization_ideal(m);
    end
    if is_zero(I)
      return (base_ring(X),X) 
    elseif is_one(I) #this means the matrix is not realizable by self-projecting points
      return (base_ring(I), nothing)
    else
      R = base_ring(X)
      QR, phi = quo(R,I)
      return (QR, phi.(X))
    end
  end
end


@doc raw"""
  selfprojecting_realization_space(m::Matroid; B::Union{GroundsetType,Nothing}=nothing; check::Bool = true)

  Function to compute the selfprojecting realization space of a selfprojecting matroid. 
  A basis B can be given that will correspond to the identity matrix in the realization. If nothing is given, a choice will be made.   
  The kwarg check states whether the matroid will be checked for self-projectivity. The default is true.

  !!! warning "This function uses the computation of the selfprojecting_realization_ideal. Therefore, it is slow except for small matroids."
# Examples
```jldoctest
julia> m = matroid_from_nonbases([[1,2,3],[4,5,6]],6)
Matroid of rank 3 on 6 elements

julia> selfprojecting_realization_space(m)
The selfprojecting realization space is
  [1   0   1   0    1    1]
  [0   1   1   0   x1   x1]
  [0   0   0   1    1   x2]
in the multivariate polynomial ring in 2 variables over QQ
avoiding the zero loci of the polynomials
RingElem[x2, -x1, -x2 + 1, -x1 + 1]

julia> m = uniform_matroid(3,6)
Matroid of rank 3 on 6 elements

julia> selfprojecting_realization_space(m,B=[4,5,6])
The selfprojecting realization space is
  [1    1    1   1   0   0]
  [1   x1   x3   0   1   0]
  [1   x2   x4   0   0   1]
in the multivariate polynomial ring in 4 variables over QQ
within the vanishing set of the ideal
Ideal (x1*x2*x3 - x1*x2*x4 - x1*x3*x4 + x1*x4 + x2*x3*x4 - x2*x3)
avoiding the zero loci of the polynomials
RingElem[x1*x4 - x1 - x2*x3 + x2 + x3 - x4, -x1 + x2, -x2 + 1, x1 - 1, -x3 + x4, -x4 + 1, x3 - 1, x1*x4 - x2*x3, x2 - x4, -x1 + x3, x2, -x1, x4, -x3]
```
"""
function selfprojecting_realization_space(m::Matroid;
  B::Union{GroundsetType,Nothing}=nothing, check::Bool = true)::MatroidRealizationSpaceSelfProjecting
  if check
    @req is_selfprojecting(m) "The given matroid is not self-projecting."
  end
  @req length(loops(m))==0 "This method is currently not implemented for matroids with loops." # as soon as realization_space is implemented for matroids with loops, this requirement can fall.
  RS = realization_space(m,char=0,simplify = true, saturate = true,ground_ring = QQ)
  R = ambient_ring(RS)
  if !is_realizable(m,char=0)
    RS = MatroidRealizationSpaceSelfProjecting(defining_ideal(RS), inequations(RS), R, nothing, 0, nothing, QQ)
    return RS
  end
  I = selfprojecting_realization_ideal(m);
  n = length(matroid_groundset(m))
  goodM = isomorphic_matroid(m, [i for i in 1:n])
  Bs = bases(goodM)
  if !isnothing(B)
    goodB = sort!(Int.([m.gs2num[j] for j in B]))
  else
    goodB = find_good_basis_heuristically(goodM)
  end
  RR, M = selfprojecting_realization_matrix(goodM, goodB, I = I) 
  if M == nothing 
    Ineqs = inequations(RS);
  else
    Ineqs = basis_minors(M,bases(m));
  end
  return MatroidRealizationSpaceSelfProjecting(I, Ineqs, R, M, 0, nothing, QQ)
end


@doc raw"""
  dimension(MRS::MatroidRealizationSpaceSelfProjecting)::Int

  Function to compute the dimension of the selfprojecting realization space of a selfprojecting matroid. 

# Examples
```jldoctest
julia> dimension(selfprojecting_realization_space(uniform_matroid(3,6)))
3
```
"""

function dimension(MRS::MatroidRealizationSpaceSelfProjecting)::Union{Int,NegInf}
  if iszero(defining_ideal(MRS))
    return length(gens(base_ring(defining_ideal(MRS))))
  else
    return dim(defining_ideal(MRS))
  end
end

@doc raw"""
  dimension(MRS::MatroidRealizationSpace)::Int

  Function to compute the dimension of the realization space of a matroid. 

# Examples
```jldoctest
julia> dimension(realization_space(uniform_matroid(3,6)))
4
```
"""
function dimension(MRS::MatroidRealizationSpace)::Union{Int,NegInf}
  if iszero(defining_ideal(MRS))
    return length(gens(base_ring(defining_ideal(MRS))))
  else
    return dim(defining_ideal(MRS))
  end
end
