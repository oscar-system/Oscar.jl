#Fix the following: KamelCase for struct and snake_case for functions! #Other naming options: SelfProjectingMRS or MatroidRealizationSpaceSelfProjecting
@attributes mutable struct MatroidRealizationSpace_SelfProj{BaseRingType, RingType} <: AbsAffineScheme{BaseRingType, RingType}
  defining_ideal::Union{Ideal,NumFieldOrderIdeal}
  inequations::Vector{RingElem}
  ambient_ring::Ring
  selfproj_realization_matrix::Union{MatElem,Nothing}
  char::Union{Int,Nothing}
  q::Union{Int,Nothing}
  ground_ring::Ring
  one_realization::Bool

# Do I need underlying_scheme for my struct?  
  # Fields for caching
  underlying_scheme::AbsAffineScheme{BaseRingType, RingType}

  function MatroidRealizationSpace_SelfProj(
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

function underlying_scheme(RS::MatroidRealizationSpace_SelfProj{BRT, RT}) where {BRT<:Ring, RT<:MPolyQuoLocRing}
  isdefined(RS, :underlying_scheme) && return RS.underlying_scheme::AffineScheme{BRT, RT}

  P = ambient_ring(RS)::MPolyRing
  I = defining_ideal(RS)::MPolyIdeal
  U = MPolyPowersOfElement(P, P.(inequations(RS)))::MPolyPowersOfElement
  RS.underlying_scheme = spec(P, I, U)
  return RS.underlying_scheme::AffineScheme{BRT, RT}
end


function Base.show(io::IO, ::MIME"text/plain", RS::MatroidRealizationSpace_SelfProj)
  if has_attribute(RS, :is_realizable) && !is_realizable(RS)
    if RS.char === nothing && RS.q === nothing
      print(io, "The matroid does not have a self-projecting realization.")
    else
      print(io, "The matroid does not have a self-projecting realization over the specified field or characteristic.")
    end
  else
    io = Oscar.pretty(io)
    if RS.one_realization
      println(io, "One selfprojecting realization is given by")
    elseif has_attribute(RS, :is_realizable) && is_realizable(RS)
      println(io, "The selfprojecting ealizations are parametrized by")
    elseif !has_attribute(RS, :is_realizable)
      println(io, "The selfprojecting realization space is")
    end
    print(io, Oscar.Indent())
    show(io, MIME("text/plain"), RS.selfproj_realization_matrix)
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
    defining_ideal(RS::MatroidRealizationSpace_SelfProj)

The ideal of the matroid realization space `RS`.
"""
defining_ideal(RS::MatroidRealizationSpace_SelfProj) = RS.defining_ideal

@doc raw"""
    inequations(RS::MatroidRealizationSpace_SelfProj)

Generators of the localizing semigroup of `RS`.
These are the polynomials that need to be nonzero in any realization.
"""
inequations(RS::MatroidRealizationSpace_SelfProj) = RS.inequations

@doc raw"""
    ambient_ring(RS::MatroidRealizationSpace_SelfProj)

The polynomial ring containing the ideal `defining_ideal_sp(RS)` and the polynomials in `inequations(RS)`.
"""
ambient_ring(RS::MatroidRealizationSpace_SelfProj) = RS.ambient_ring



# #Can we cahnge the following to a function for selfrpojecting_realizations?
# @doc raw"""
#     is_realizable(M; char::Union{Int,Nothing}=nothing, q::Union{Int,Nothing}=nothing)

# * If char = nothing, then this function determines whether the matroid is realizable over some field.

# * If `char == 0`, then this function determines whether the matroid is realizable over some field of
#     characteristic 0.

# * If char = p is prime, this function determines whether the matroid is realizable
#     over the finite field ``GF(p)``.

# * If `char == p` and `q` is a power of `p`, this function determines whether the matroid is realizable over the
#     finite field ``GF(q)``.
# """
# function is_realizable_sp(
#   M::Matroid; char::Union{Int,Nothing}=nothing, q::Union{Int,Nothing}=nothing
# )::Bool
#   RS = realization_space(M; char=char, q=q)
#   return is_realizable_sp(RS)
# end

# @attr Bool function is_realizable_sp(RS::MatroidRealizationSpace_SelfProj)
#   if !(RS.ambient_ring_sp isa MPolyRing)
#     return true
#   end
#   for p in minimal_primes(RS.defining_ideal_sp)
#     component_non_trivial = all(!in(p), RS.inequations_sp)
#     if component_non_trivial
#       return true
#     end
#   end
#   return false
# end


@doc raw"""
    selfproj_realization_matrix(RS::MatroidRealizationSpace_SelfProj)

A matrix with entries in ambient_ring_sp(RS) whose columns, when filled in with values satisfying equalities
from `defining_ideal(RS)` and inequations from `inequations(RS)`, form a self-projecting realization for the matroid.
"""
selfproj_realization_matrix(RS::MatroidRealizationSpace_SelfProj) = RS.selfproj_realization_matrix

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
    is_selfprojecting(mat::Matroid)

Returns a boolean which states whether the given matroid satisfies the property to be self-projecting.
# Examples
```jldoctest
julia> M = fano_matroid();
```
"""

function selfproj_realization_ideal(m::Matroid;saturate::Bool = false)::Ideal
    @warn "This function is slow for...?"
    if !is_selfprojecting(m) 
      error("The given matroid is not self-projecting.") #Is it too costly to have this check?! How about an option to check, which can be turned off?
    end
    if !is_realizable(m)
      #What to do? Error? or "empty" What is the default for MRS?
    end
    RS = realization_space(m,char =0,simplify = true, saturate = true); 
    R = ambient_ring(RS); #I can't call on the variables of R.
    I = defining_ideal(RS);
    n = length(matroid_groundset(m))
    k = rank(m)
    RR, x = polynomial_ring(QQ, :x=> 1:length(gens(R)))
    F = hom(R, RR,gens(RR)) #use this to move ideals into RR
    G = hom(RR, R, gens(R)) #use this to move them back again into the ambient ring of RS
    S, l = polynomial_ring(RR, :l => 1:n)
    L = ideal(S, prod(l[i] for i in 1:n))
    MRS = realization_matrix(RS)
    M = matrix(RR,nrows(MRS),ncols(MRS),[F(MRS[i,j]) for i in 1:nrows(MRS) for j in 1:ncols(MRS)])
    DL = diagonal_matrix(S,[l[i] for i in 1:n])
    V = M*DL*transpose(M)
    IS = ideal(S,[V[i,j] for i in 1:k for j in 1:k]) + ideal(S,[F(gens(I)[i]) for i in 1:length(gens(I))])
    J = stepwise_saturation(IS,l) 
    E = eliminate(J,l) 
    if saturate 
      E = stepwise_saturation(E,inequations(RS))
    end
    if is_zero(E)
      return ideal(ambient_ring(RS),[0]);
    elseif is_one(E)
      return ideal(ambient_ring(RS),[1]);
    else
      return ideal(R,G.(gens(E)))
    end
end
@doc raw"""
    selfproj_realization_ideal(m::Matroid)

    Function to compute the defining_ideal of a selfprojecting realization space
"""

#this function computes the minors of a given matrix for the bases in a given list.
function basisminors(M::MatElem, Bases::Vector{Vector{Int}})::Vector
    R = base_ring(M);
    return [R(det(M[1:nrows(M),Bases[i]])) for i in 1:length(Bases)]
end

###################
#B are the given columns that will be the identity matrix
# it might be easier for me to just take the standard realization matrix, and then use a homomorphism into the quotient ring by the selfproj_realization_ideal to simplify (check how the simplify functions work!)
function selfproj_realization_matrix(M::Matroid, B::Vector{Int}, F::Ring)
  #include a check that M is realizable & selfproj
  if !is_selfprojecting(M) 
    error("The given matroid is not self-projecting.")
  end
  if !(is_realizable(M))
    return nothing
  end
  RSM = realization_space_matrix(M,B,F)
  X = RSM[2];
  I = selfproj_realization_ideal(M);
  if is_zero(I)
    return X
  end
  if is_one(I) #this means the matrix is not realizable by self-projecting points
    return nothing
  end
  #need to test this!
  R = RSM[1] 
  xs = gens(R)
  cR = coefficient_ring(R)
  nr, nc = size(X)
  QR, phi = quo(R,I)
  return phi.(X)
end


function selfprojecting_realization_space(m::Matroid;
  B::Union{GroundsetType,Nothing}=nothing)::MatroidRealizationSpace_SelfProj
  if !is_selfprojecting(m) 
    error("The given matroid is not self-projecting.")
  end
  RS = realization_space(m,char=0,simplify = true, saturate = true)
  R = ambient_ring(RS)
  I = selfproj_realization_ideal(m);
  n = length(matroid_groundset(m))
  goodM = isomorphic_matroid(m, [i for i in 1:n])
  Bs = bases(goodM)
  if !isnothing(B)
    goodB = sort!(Int.([M.gs2num[j] for j in B]))
  else
    goodB = find_good_basis_heuristically(goodM)
  end
  M = selfproj_realization_matrix(goodM, goodB, RS.ground_ring) #do I want selfproj_realization_matrix to returna tuple of ring and matrix? Currently just the matrix, since this is what i need most urgently
  Ineqs = basisminors(M,bases(m));
  ##include a check for realizable, how does the set_attribute!(MRS, :is_realizable, :false) even work for me?

  #M = selfproj_realization_matrix(m,b,R)

  ##need to find reference basis b for the identity matrix in the realizatio matrix
  ##R is ambient_ring(RS) where RS is the standard realization space
  ##ground_ring hardcoded as QQ
  ##do I need boo the boolean for one_realization_sp?
  return MatroidRealizationSpace_SelfProj(I, Ineqs, R, M, 0, nothing, QQ)
  
end




#the 2 functions below need to be tested 
function dimension(MRS::MatroidRealizationSpace_SelfProj)::Int
  return dim(defining_ideal(MRS))
end

function dimension(MRS::MatroidRealizationSpace)::Int
  return dim(defining_ideal(MRS))
end

# Exports
# export MatroidRealizationSpace_SelfProj
# export inequations
# export defining_ideal
# export underlying_scheme
# export selfproj_realization_ideal
# export dimension
