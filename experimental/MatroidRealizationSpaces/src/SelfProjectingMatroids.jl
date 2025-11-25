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

#I need currently have no analogue for is_realizable and one_realization for MatroidRealizationSpaceSelfProjecting
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
These are the polynomials that need to be nonzero in any realization.
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
# Examples
```jldoctest
julia> m = matroid_from_nonbases([[1,2,3],[4,5,6]],6)
Matroid of rank 3 on 6 elements

julia> selfprojecting_realization_ideal(m)
┌ Warning: This function is slow except for small matroids!
└ @ Oscar ~/Documents/Oscar.jl/experimental/MatroidRealizationSpaces/src/SelfProjectingMatroids.jl:169
Ideal generated by
  0

julia> m = uniform_matroid(3,6)
Matroid of rank 3 on 6 elements

julia> selfprojecting_realization_ideal(m)
┌ Warning: This function is slow except for small matroids!
└ @ Oscar ~/Documents/Oscar.jl/experimental/MatroidRealizationSpaces/src/SelfProjectingMatroids.jl:169
Ideal generated by
  x1*x2*x3 - x1*x2*x4 - x1*x3*x4 + x1*x4 + x2*x3*x4 - x2*x3
  ```
"""
#this function is not properly tested, since it did not terminate for intersting examples.
function selfprojecting_realization_ideal(m::Matroid;saturate::Bool = false, check::Bool = true)::Ideal
    @warn "This function is slow except for small matroids!"
    if check 
      @req is_self_projecting(m) "The given matroid is not self projecting" 
    end
    if !is_realizable(m,char = 0)
      return defining_ideal(realization_space(m,char =0, QQ))
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
    IS = ideal(RR,[V[i,j] for i in 1:k for j in 1:k]) + ideal(RR,[F(gens(I)[i]) for i in 1:length(gens(I))])
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


#this function computes the minors of a given matrix for the bases in a given list.
#this list is not reduced, i.e. there are no repetitions but the polynomials are not reduced, i.e. the list might contain x+1 and y-2 and their product
function sort_by_degree(polys::Vector{<:RingElem})::Vector{<:RingElem}
    deg(p) = hasmethod(total_degree, Tuple{typeof(p)}) ? total_degree(p) : 0
    return sort(polys, by=deg)
end

function basisminors(M::MatElem, Bases::Vector{Vector{Int}})::Vector{<:RingElem}
    R = base_ring(M);
    candidates = sort_by_degree([R(det(M[1:nrows(M),Bases[i]])) for i in 1:length(Bases)])
    multiplicativeSet= nothing
    ineqs = [R(0)]
    if R isa MPolyQuoRing #in this case I cannot make the multiplicativeSet using the powers_of_element so we currently use the lesser choice: here there might be double entries or products of polynomials in the list - is there a better solution for this?
      for i in 1:length(candidates)
        @req !iszero(candidates[i]) "a basis has vanishing minor"
        if isone(candidates[i]) || isone(-candidates[i]) 
        elseif !(candidates[i] in ineqs[2:length(ineqs)])&& !(-candidates[i] in ineqs[2:length(ineqs)])
          push!(ineqs,R(candidates[i]))
        end
      end
      return ineqs[2:length(ineqs)]
    end
    for i in 1:length(candidates)
      if isone(candidates[i]) || isone(-candidates[i]) 
      elseif iszero(candidates[i])
        error("a basis has vanishing minor")
      else
        if isnothing(multiplicativeSet) && !(candidates[i] in ineqs[2:length(ineqs)]) && !(-candidates[i] in ineqs[2:length(ineqs)])
          multiplicativeSet = powers_of_element(candidates[i])
          push!(ineqs,R(candidates[i]))
        elseif !(candidates[i] in ineqs[2:length(ineqs)])&& !(-candidates[i] in ineqs[2:length(ineqs)]) && !(candidates[i] in multiplicativeSet)  && !(-candidates[i] in multiplicativeSet)# !(candidates[i] in ineqs[2:length(ineqs)])&& !(-candidates[i] in ineqs[2:length(ineqs)])
          #one possibility is to take out all products of existing ineqs
          multiplicativeSet = product(multiplicativeSet,powers_of_element(candidates[i]));
          push!(ineqs,R(candidates[i]))
        end
      end
    end
    return ineqs[2:length(ineqs)]
end

###################
#Bas are the given columns that will be the identity matrix
# it might be easier to just take the standard realization matrix, and then use a homomorphism into the quotient ring by the selfprojecting_realization_ideal to simplify (check how the simplify functions work!)
#this function is not properly tested, since it did not terminate for intersting examples.
function selfprojecting_realization_matrix(M::Matroid, Bas::Vector{Int}, F::Ring)
  #include a check that M is realizable & selfproj
  if !is_selfprojecting(M) 
    error("The given matroid is not self-projecting.")
  end
  if !(is_realizable(M))
    return nothing
  end
  #RSM
  X = realization_matrix(realization_space(M,B=Bas,ground_ring=QQ)) #this matrix should be simplified
  #X = RSM[2];
  I = selfprojecting_realization_ideal(M);
  if is_zero(I)
    return X #In this case S = R, I need to use the simplified X to make the output nice
  elseif is_one(I) #this means the matrix is not realizable by self-projecting points
    return nothing
  else
  R = base_ring(X)
 # xs = gens(R)
 # cR = coefficient_ring(R)
 # nr, nc = size(X)
  QR, phi = quo(R,I)
  return phi.(X)
  end
end

#function to compute the self-projecting realization space
function selfprojecting_realization_space(m::Matroid;
  B::Union{GroundsetType,Nothing}=nothing)::MatroidRealizationSpaceSelfProjecting
  if !is_selfprojecting(m) 
    error("The given matroid is not self-projecting.")
  end
  if length(loops(m))>0
    error("This method is currently not implemented matroids with loops.")
  end
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
    goodB = sort!(Int.([M.gs2num[j] for j in B]))
  else
    goodB = find_good_basis_heuristically(goodM)
  end
  M = selfprojecting_realization_matrix(goodM, goodB, RS.ground_ring) #does not return a tuple of ring and matrix like it does for realization_space #inorder to avoid calling selfprojecting_realization_ideal twice, one could modify the function to give it the already computed ideal, maybe optionally?
  if M == nothing 
    Ineqs = inequations(RS);
  else
    Ineqs = basisminors(M,bases(m));
  end
  return MatroidRealizationSpaceSelfProjecting(I, Ineqs, R, M, 0, nothing, QQ)
end

#do we need this function to make something faster? or should it just be included in case someone wants to test something?
function veronese2(M::MatElem)
    R = base_ring(M)
    newCols = Vector{Vector{elem_type(R)}}()
    for j in 1:ncols(M)
        v = M[:, j]
        col = elem_type(R)[]
        for i in 1:nrows(M)
            for k in i:nrows(M)
                push!(col, v[i] * v[k])
            end
        end
        push!(newCols, col)
    end
    return transpose(matrix(R, hcat(newCols...)))
end

#the 2 functions below need to be tested, they do not seem to give the correct answer!
# function dimension(MRS::MatroidRealizationSpaceSelfProjecting)::Int
#   return dim(defining_ideal(MRS))
# end

# function dimension(MRS::MatroidRealizationSpace)::Int
#   return dim(defining_ideal(MRS))
# end
