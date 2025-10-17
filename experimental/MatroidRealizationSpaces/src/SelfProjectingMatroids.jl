#Fix the following: KamelCase for struct and snake_case for functions!
@attributes mutable struct MatroidRealizationSpace_SelfProj{BaseRingType, RingType} <: AbsAffineScheme{BaseRingType, RingType}
  defining_ideal::Union{Ideal,NumFieldOrderIdeal}
  inequations::Vector{RingElem}
  ambient_ring::Ring
  realization_matrix::Union{MatElem,Nothing}
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
    show(io, MIME("text/plain"), RS.realization_matrix)
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



##Can we cahnge the following to a function for selfrpojecting_realizations?
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


# function underlying_scheme(RS::MatroidRealizationSpace_SelfProj)
#   error("method not implemented for realization spaces of type $(typeof(RS))")
# end
# @doc raw"""
#     realization_matrix_sp(RS::MatroidRealizationSpace_SelfProj)

# A matrix with entries in ambient_ring_sp(RS) whose columns, when filled in with values satisfying equalities
# from `defining_ideal_sp(RS)` and inequations_sp from `inequations_sp(RS)`, form a realization for the matroid.
# """
# realization_matrix_sp(RS::MatroidRealizationSpace_SelfProj) = RS.realization_matrix_sp

#A function to test if a matroid is self-projecting
function is_selfprojecting(mat::Matroid)::Bool
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
julia> M = fano_matroid();
```
"""

#how to fix the error with returning the ideal? I do want the ambient_rings of MRS and MRS_SP to be compatible!
function selfproj_realization_ideal(m::Matroid)::MatroidRealizationSpace_SelfProj
    print("Warning: This function is slow for...?")
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
    J = saturation(IS,L)
    E = eliminate(J,l) 
    if gens(E)[1] == 0
      return ideal(ambient_ring(RS),[0]);
      # Cannot convert an object of type MPolyIdeal{QQMPolyRingElem} to an object of type MatroidRealizationSpace_SelfProj
    else
      return ideal(R,[G(gens(E)[i]) for i in 1:length(gens(E))])
    end
end
@doc raw"""
    selfproj_realization_ideal(m::Matroid)

    Function to compute the defining_ideal of a selfprojecting realization space
"""

function selfprojecting_realization_space(m::Matroid)::MatroidRealizationSpace_SelfProj
  if !is_selfprojecting(m) 
    error("The given matroid is not self-projecting.")
  end
  RS = realization_space(m,char=0,simplify = true, saturate = true)
  R = ambient_ring(RS)
  #need to find reference basis b for basisminors
  #R is ambient_ring(RS) where RS is the standard realization space
  #ground_ring hardcoded as QQ
  #do I need boo the boolean for one_realization_sp?
  return MatroidRealizationSpace_SelfProj(selfproj_realization_ideal(m), basisminors(m,b), R, selfproj_realization_matrix(m), 0, nothing, QQ, boo)
end


###################
# function realization_space_matrix(M::Matroid, B::Vector{Int}, F::Ring)
#   # prepare the combinatorial data

#   circs = fundamental_circuits_of_basis(M, B)

#   nonIdCols = setdiff(matroid_groundset(M), B)
#   circs = [setdiff(c, nonIdCols) for c in circs]

#   rk = rank(M)
#   n = length(M)

#   # we start by computing the number of variables:
#   numVars = 0
#   unUsedRowsForOnes = collect(2:rk)
#   for col in 1:(n - rk), row in 1:rk
#     circ = circs[col]
#     if !(B[row] == minimum(circ)) && B[row] in circ
#       if row in unUsedRowsForOnes
#         unUsedRowsForOnes = setdiff(unUsedRowsForOnes, [row])
#       else
#         numVars += 1
#       end
#     end
#   end

#   if numVars > 0
#     R, x = polynomial_ring(F, numVars)
#   else
#     R = F
#     x = Vector{MPolyRingElem}()
#   end

#   unUsedRowsForOnes = collect(2:rk)

#   # create the matrix and fill it with entries

#   mat = zero_matrix(R, rk, n)

#   for i in 1:rk
#     mat[i, B[i]] = R(1)
#   end

#   varCounter = 1

#   for col in 1:(n - rk), row in 1:rk
#     circ = circs[col]
#     c = nonIdCols[col]

#     if B[row] == minimum(circ)
#       mat[row, c] = R(1)
#     elseif B[row] in circ
#       if row in unUsedRowsForOnes
#         mat[row, c] = R(1)
#         unUsedRowsForOnes = setdiff(unUsedRowsForOnes, [row])
#       else
#         mat[row, c] = x[varCounter]
#         varCounter = varCounter + 1
#       end
#     else
#       mat[row, c] = R(0)
#     end
#   end
#   return (R, mat)
# end


# @doc raw"""
#     realization_space(
#       M::Matroid;
#       B::Union{GroundsetType,Nothing}=nothing,
#       saturate::Bool=false,
#       simplify::Bool=true,
#       char::Union{Int,Nothing}=nothing,
#       q::Union{Int,Nothing}=nothing,
#       ground_ring::Ring=ZZ
#     )::MatroidRealizationSpace

# This function returns the data for the coordinate ring of the matroid realization space of the matroid `M`
# as a `MatroidRealizationSpace`. This function has several optional parameters.

# * `B` is a basis of M that specifies which columns of `realization_matrix_sp(M)` form an identity matrix.
#   The default is `nothing`, in which case the basis is chosen for you.

# * `saturate` determines whether `defining_ideal_sp(M)` should be saturated with respect to the semigroup
#   generated by `inequations_sp(M)`. The default is `false`. The saturation can be rather slow for large instances.

# * `simplify` determines whether a reduced realization space is returned which means that the equations
#   are used to eliminate variables as far as possible. The default is `true`.

# * `char` specifies the characteristic of the coefficient ring. The returned realization space is then the space of all
#   realizations over fields of characteristic `char`. The default is `nothing`.

# * `q` is an integer and assumed to be a prime power `q=p^k`. The returned realization space is then the space of all
#   realizations over the field ``GF(p^k)``. The default is `nothing`.

# * `ground_ring` is a ring and specifies the ground_ring over which one wants to consider the realization space,
#   e.g. `QQ` or `GF(p)`. The groud_ring `ZZ` means that we compute the space of realizations over all fields.
#   The default is `ZZ`.

# # Examples
# ```jldoctest
# julia> M = fano_matroid();

# julia> RS = realization_space(M)
# The realization space is
#   [0   1   1   1   1   0   0]
#   [1   0   1   1   0   1   0]
#   [1   0   1   0   1   0   1]
# in the integer ring
# within the vanishing set of the ideal
# 2ZZ

# julia> realization_space(non_fano_matroid())
# The realization space is
#   [1   1   0   0   1   1   0]
#   [0   1   1   1   1   0   0]
#   [0   1   1   0   0   1   1]
# in the integer ring
# avoiding the zero loci of the polynomials
# RingElem[2]

# julia> realization_space(pappus_matroid(), char=0)
# The realization space is
#   [1   0   1   0   x2   x2                 x2^2    1    0]
#   [0   1   1   0    1    1   -x1*x2 + x1 + x2^2    1    1]
#   [0   0   0   1   x2   x1                x1*x2   x1   x2]
# in the multivariate polynomial ring in 2 variables over QQ
# avoiding the zero loci of the polynomials
# RingElem[x1 - x2, x2, x1, x2 - 1, x1 + x2^2 - x2, x1 - 1, x1*x2 - x1 - x2^2]

# julia> realization_space(uniform_matroid(3,6))
# The realization space is
#   [1   0   0   1    1    1]
#   [0   1   0   1   x1   x3]
#   [0   0   1   1   x2   x4]
# in the multivariate polynomial ring in 4 variables over ZZ
# avoiding the zero loci of the polynomials
# RingElem[x1*x4 - x2*x3, x2 - x4, x1 - x3, x1*x4 - x1 - x2*x3 + x2 + x3 - x4, x3 - x4, x4 - 1, x3 - 1, x3, x4, x1 - x2, x2 - 1, x1 - 1, x1, x2]
# ```
# """
# function realization_space(
#   M::Matroid;
#   B::Union{GroundsetType,Nothing}=nothing,
#   saturate::Bool=false,
#   simplify::Bool=true,
#   char::Union{Int,Nothing}=nothing,
#   q::Union{Int,Nothing}=nothing,
#   ground_ring::Ring=ZZ
# )::MatroidRealizationSpace
#   if char != nothing && !is_prime(char) && char != 0
#     error("The characteristic has to be 0 or a prime number.")
#   end

#   #Construct the base ring as F_p if q=p^k
#   if q != nothing
#     isprimepower, k, p = is_prime_power_with_data(q)
#     if !isprimepower
#       error("The given q has to be a prime power.")
#     end
#     if char != nothing && char != p
#       error("The given characteristic doesn't match q.")
#     else
#       char = p
#     end
#   end

#   if char == 0
#     ground_ring = QQ
#   elseif char != nothing
#     ground_ring = GF(char)
#   end

#   rk = rank(M)
#   n = length(M)

#   goodM = isomorphic_matroid(M, [i for i in 1:n])

#   Bs = bases(goodM)

#   if !isnothing(B)
#     goodB = sort!(Int.([M.gs2num[j] for j in B]))
#   else
#     goodB = find_good_basis_heuristically(goodM)
#   end
#   polyR, mat = realization_space_matrix(goodM, goodB, ground_ring)

#   eqs = Vector{RingElem}()
#   ineqs = Vector{RingElem}()

#   #need to catch the corner-case if there are no variables at all
#   if !(polyR isa MPolyRing)
#     RS = MatroidRealizationSpace(ideal(polyR, [0]), ineqs, polyR, mat, char, q, ground_ring)
#     set_attribute!(RS, :is_realizable_sp, :true)
#     return RS
#   end

#   for col in subsets(Vector(1:n), rk)
#     col_det = det(mat[:, col])

#     if total_degree(col_det) <= 0
#       if col_det != 0 && col in Bs
#         if is_unit(col_det)
#           continue
#         end
#       elseif col_det != 0 # and col is not a basis
#         # determinant nonzero but set not a basis
#         push!(eqs, col_det)
#       elseif col in Bs
#         push!(ineqs, col_det)
#         #determinant zero but set is a basis, i.e. M is not realizable
#         RS = MatroidRealizationSpace(ideal(polyR, eqs), ineqs, polyR, nothing, char, q, ground_ring)
#         set_attribute!(RS, :is_realizable_sp, :false)
#         return RS
#       else
#         continue
#       end
#     end

#     if col in Bs
#       push!(ineqs, col_det)
#     else
#       push!(eqs, col_det)
#     end
#   end

#   def_ideal = ideal(polyR, eqs)
#   def_ideal = ideal(groebner_basis(def_ideal))
#   if isone(def_ideal)
#     RS = MatroidRealizationSpace(def_ideal, ineqs, polyR, nothing, char, q, ground_ring)
#     set_attribute!(RS, :is_realizable_sp, :false)
#     return RS
#   end

#   ineqs = gens_2_prime_divisors(ineqs)

#   RS = MatroidRealizationSpace(def_ideal, ineqs, polyR, mat, char, q, ground_ring)

#   if simplify
#     RS = reduce_realization_space(RS)
#   end

#   if q != nothing && RS.ambient_ring_sp isa MPolyRing
#     I = RS.defining_ideal_sp
#     R = RS.ambient_ring_sp
#     eqs = Vector{RingElem}()
#     for x in gens(R)
#       push!(eqs, x^q - x)
#     end
#     I = I + ideal(R, eqs)
#     RS.defining_ideal_sp = I
#     if isone(RS.defining_ideal_sp)
#       set_attribute!(RS, :is_realizable_sp, :false)
#       return RS
#     end
#   end

#   if saturate
#     RS.defining_ideal_sp = stepwise_saturation(RS.defining_ideal_sp, RS.inequations_sp)
#     if isone(RS.defining_ideal_sp)
#       set_attribute!(RS, :is_realizable_sp, :false)
#       return RS
#     end
#   end

#   return RS
# end

# # A heuristic function that tries to find a sensible basis for the moduli space computation for which the defining ideal is not too complicated
# function find_good_basis_heuristically(M::Matroid)
#   bs = bases(M)
#   cs = circuits(M)
#   min_num_vars = length(cs) * rank(M)
#   min_basis = bs[1]
#   for b in bs
#     current_num_vars = 0
#     for c in cs, e in c
#       if e in b
#         current_num_vars += 1
#       end
#     end
#     if current_num_vars < min_num_vars
#       min_num_vars = current_num_vars
#       min_basis = b
#     end
#   end
#   return min_basis
# end

# # Return the prime divisors of f.
# function poly_2_prime_divisors(f::RingElem)
#   return map(first, factor(f))
# end

# # Return the unique prime divisors of the elements of Sgen, again no exponents.
# function gens_2_prime_divisors(Sgens::Vector{<:RingElem})
#   return unique!(vcat([poly_2_prime_divisors(f) for f in Sgens]...))
# end

# @doc raw"""
#     realization(M::Matroid; B::Union{GroundsetType,Nothing} = nothing,
#       saturate::Bool=false,
#       char::Union{Int,Nothing}=nothing, q::Union{Int,Nothing}=nothing
#     )::MatroidRealizationSpace

# This function tries to find one realization in the matroid realization space of
# the matroid `M`. The output is again a `MatroidRealizationSpace`.

# If the matroid is only realizable over an extension of the prime field the
# extension field is specified as a splitting field of an irreducible polynomial.
# Every root of this polynomial gives an equivalent realization of the matroid.

# This function has several optional parameters. Note that one must input either
# the characteristic or a specific field of definition for the realization.

# * `B` is a basis of M that specifies which columns of `realization_matrix_sp(M)`
#   form the identity matrix. The default is `nothing`, in which case the basis is
#   chosen for you.

# * `char` specifies the characteristic of the coefficient ring, and is used to determine if the matroid
#   is realizable over a field of this characteristic. The default is `nothing`.

# * `q` is an integer, and when char = p, this input is used to determine whether the matroid
#   is realizable over the finite field ``GF(p^{q})``. The default is `nothing`.

# * `reduce` determines whether a reduced realization space is returned which means that the equations
#   are used to eliminate variables as far as possible. The default is `true`.

# * `saturate` determines whether `defining_ideal_sp(M)` should be saturated with respect to the semigroup
#   generated by `inequations_sp(M)`. The default is `false`. This can be rather slow for large instances.

# # Examples
# ```jldoctest
# julia> realization(pappus_matroid(), char=0)
# One realization is given by
#   [1   0   1   0   2   2   4   1   0]
#   [0   1   1   0   1   1   1   1   1]
#   [0   0   0   1   2   3   6   3   2]
# in the rational field

# julia> realization(pappus_matroid(), q=4)
# One realization is given by
#   [1   0   1   0   x1 + 1   x1 + 1   x1    1        0]
#   [0   1   1   0        1        1    1    1        1]
#   [0   0   0   1   x1 + 1       x1    1   x1   x1 + 1]
# in the multivariate polynomial ring in 1 variable over GF(2)
# within the vanishing set of the ideal
# Ideal (x1^2 + x1 + 1)

# julia> realization(uniform_matroid(3,6), char=5)
# One realization is given by
#   [1   0   0   1   1   1]
#   [0   1   0   1   4   3]
#   [0   0   1   1   3   2]
# in the prime field of characteristic 5
# ```
# """
# function realization(
#   M::Matroid;
#   B::Union{GroundsetType,Nothing}=nothing,
#   saturate::Bool=false,
#   simplify::Bool=true,
#   char::Union{Int,Nothing}=nothing,
#   q::Union{Int,Nothing}=nothing,
# )
#   RS = realization_space(M; B=B, saturate=saturate, simplify=simplify, char=char, q=q)

#   return realization(RS)
# end

# @doc raw"""
#     realization(RS::MatroidRealizationSpace)

# This function tries to find one realization in the matroid realization `RS`.
# The output is again a `MatroidRealizationSpace`.
# """
# function realization(RS::MatroidRealizationSpace)
#   # If the matroid is not realizable we stop
#   if RS.char == nothing && RS.q == nothing
#     error("A field or characteristic must be specified")
#   end

#   if !is_realizable_sp(RS)
#     return RS
#   end

#   # If the ambient ring is not a polynomial ring we can't reduce and we stop
#   R = RS.ambient_ring_sp

#   !(R isa MPolyRing) && return RS
#   Inew = RS.defining_ideal_sp
#   eqs = copy(gens(Inew))

#   if dim(Inew) == 0
#     for p in minimal_primes(Inew)
#       if !any(in(p), RS.inequations_sp)
#         Inew = p
#         break
#       end
#     end
#     RSnew = MatroidRealizationSpace(Inew, Vector{RingElem}(), R, RS.realization_matrix_sp, RS.char, RS.q, RS.ground_ring)
#     RSnew = reduce_realization_space(RSnew)
#     RSnew.one_realization_sp = true
#     return RSnew
#   end

#   d = min(dim(Inew), nvars(R))
#   ineqsnew = RS.inequations_sp

#   counter = 0
#   base = 7
#   if RS.char != nothing && RS.char > 0
#     base = RS.char
#   end
#   upperbound = min(base^d, 10^3)
#   while counter < upperbound
#     eqsnew = copy(eqs)
#     values = digits(counter; base=base, pad=d)
#     counter += 1

#     for i in 1:d
#       push!(eqsnew, R[i] - values[i])
#     end
#     Inew = ideal(groebner_basis(ideal(R, eqsnew)))
#     isone(Inew) && continue
#     ineqsnew = copy(RS.inequations_sp)
#     for i in 1:length(ineqsnew)
#       ineqsnew[i] = reduce(ineqsnew[i], gens(Inew))
#     end
#     (R(0) in ineqsnew) && continue
#     Inew = stepwise_saturation(Inew, ineqsnew)
#     isone(Inew) && continue
#     break
#   end
#   counter == upperbound && d != 0 && return RS

#   RSnew = MatroidRealizationSpace(Inew, ineqsnew, R, RS.realization_matrix_sp, RS.char, RS.q, RS.ground_ring)
#   RSnew = reduce_realization_space(RSnew)
#   ineqsnew = RSnew.inequations_sp
#   if length(ineqsnew) > 0
#     Inew = RSnew.defining_ideal_sp
#     Rnew = RSnew.ambient_ring_sp
#     ineqsnew = filter(p -> !isone(ideal(groebner_basis(Inew + ideal(Rnew, p)))), ineqsnew)
#     RSnew = MatroidRealizationSpace(
#       Inew, ineqsnew, Rnew, RSnew.realization_matrix_sp, RSnew.char, RSnew.q, RSnew.ground_ring
#     )
#   end

#   RSnew.one_realization_sp = true

#   return RSnew
# end

# #####################
# # full reduction    #
# #####################

# function coefficient_v(v::RingElem, f::RingElem)
#   isone(degree(f, v)) || return (false, parent(f)(0))
#   return (true, coeff(f, [v], [1]))
# end

# function find_solution_v(
#   v::RingElem, Igens::Vector{<:RingElem}, Sgens::Vector{<:RingElem}, R::MPolyRing
# )
#   with_v_deg_1 = [g for g in Igens if isone(degree(g, v))]
#   length(with_v_deg_1) != 0 || return (false, R(0))

#   for f in with_v_deg_1
#     (a, den) = coefficient_v(v, f)
#     a || return (false, R(0))
#     fac_den = poly_2_prime_divisors(den)
#     !issubset(fac_den, Sgens) && continue
#     no_v = coeff(f, [v], [0])
#     iszero(length(no_v)) && continue
#     h = R(-1) * no_v
#     return (true, h//den)
#   end
#   return (false, R(0))
# end

# # v is replaced by t in f
# function sub_map(v::RingElem, t::RingElem, R::MPolyRing, xs::Vector{<:RingElem})
#   xs_v = map(x -> x == v ? t : x, xs)
#   return hom(R, fraction_field(R), identity, xs_v)
# end

# # replace v by t in f, only return the numerator.
# function sub_v(v::RingElem, t::RingElem, f::RingElem, R::Ring, xs::Vector{<:RingElem})
#   m = sub_map(v, t, R, xs)
#   new_f = numerator(m(f))
#   return new_f
# end

# # removes factors that are in the semigroup generated by Sgens
# function clean(f::RingElem, R::MPolyRing, Sgens::Vector{<:RingElem})
#   is_zero(f) && return f # TODO: move to the right place
#   fFactors = factor(f)
#   cleanf_arr = [k^e for (k, e) in fFactors if !(k in Sgens) || is_unit(k)]
#   length(cleanf_arr) > 0 ? prod(cleanf_arr) : unit(fFactors)
# end

# # variables in ideal
# function ideal_vars(Igens::Vector{<:RingElem})
#   return unique!(vcat([vars(gen) for gen in Igens]...))
# end

# function n_new_Sgens(
#   x::RingElem, t::RingElem, Sgens::Vector{<:RingElem}, R::Ring, xs::Vector{<:RingElem}
# )
#   preSgens = unique!([sub_v(x, t, f, R, xs) for f in Sgens])
#   if R(0) in preSgens
#     return [R(0)]
#   end
#   return gens_2_prime_divisors(preSgens)
# end

# function n_new_Igens(
#   x::RingElem,
#   t::RingElem,
#   Igens::Vector{<:RingElem},
#   Sgens::Vector{<:RingElem},
#   R::Ring,
#   xs::Vector{<:RingElem},
# )
#   preIgens = unique!([clean(sub_v(x, t, f, R, xs), R, Sgens) for f in Igens])
#   return filter(!iszero, preIgens)
# end

# function matrix_clear_den_in_col(X::Oscar.MatElem, c::Int)
#   Xc = [denominator(f) for f in X[:, c]]
#   t = lcm(Xc)
#   result = multiply_column!(X, t, c)
#   return result
# end

# function matrix_clear_den(X::Oscar.MatElem)
#   rs, cs = size(X)
#   for c in 1:cs
#     X = matrix_clear_den_in_col(X, c)
#   end
#   return X
# end

# function reduce_ideal_one_step(
#   MRS::MatroidRealizationSpace, elim::Vector{<:RingElem}, fullyReduced::Bool
# )
#   Igens = gens(MRS.defining_ideal_sp)
#   Sgens = MRS.inequations_sp
#   R = MRS.ambient_ring_sp
#   FR = fraction_field(R)
#   xs = gens(R)
#   X = MRS.realization_matrix_sp
#   nr, nc = size(X)

#   Ivars = ideal_vars(Igens)

#   t = R(0)
#   for x in Ivars
#     a, t = find_solution_v(x, Igens, Sgens, R)
#     a || continue

#     phi = sub_map(x, t, R, xs)
#     Sgens_new = n_new_Sgens(x, t, Sgens, R, xs)
#     if length(Sgens_new) == 0
#       Sgens_new = Vector{RingElem}()
#     elseif R(0) in Sgens_new
#       # 0 is in the inequations_sp, hence M is not realizable.
#       MRS.inequations_sp = Sgens_new
#       return (MRS, elim, true)
#     end

#     Igens_new = n_new_Igens(x, t, Igens, Sgens_new, R, xs)
#     push!(elim, x)

#     phiX = matrix(FR, [phi(X[i, j]) for i in 1:nr, j in 1:nc])
#     nX_FR = matrix_clear_den(phiX)
#     nX = matrix(R, [numerator(nX_FR[i, j]) for i in 1:nr, j in 1:nc])

#     GBnew = collect(groebner_basis(ideal(R, Igens_new)))

#     MRS_new = MatroidRealizationSpace(ideal(R, GBnew), Sgens_new, R, nX, MRS.char, MRS.q, MRS.ground_ring)

#     return (MRS_new, elim, fullyReduced)
#   end

#   return (MRS, elim, true)
# end

# function reduce_realization_space(
#   MRS::MatroidRealizationSpace,
#   elim::Vector{RingElem}=Vector{RingElem}(),
#   fullyReduced::Bool=false,
# )

#   #If there are no variables left, we don't reduce anything
#   if !(MRS.ambient_ring_sp isa MPolyRing)
#     return MRS
#   end

#   (MRS, elim, fullyReduced) = reduce_ideal_one_step(MRS, elim, fullyReduced)

#   !fullyReduced && return reduce_realization_space(MRS, elim, fullyReduced)

#   R = MRS.ambient_ring_sp
#   xs = gens(R)
#   cR = coefficient_ring(R)
#   X = MRS.realization_matrix_sp
#   nr, nc = size(X)
#   Igens = gens(MRS.defining_ideal_sp)
#   Sgens = MRS.inequations_sp

#   # 0 is in the inequations_sp, thus it is not realizable
#   if R(0) in Sgens
#     MRS.realization_matrix_sp = nothing
#     set_attribute!(MRS, :is_realizable_sp, :false)
#     return MRS
#   end
#   xnew_str = ["x$i" for i in 1:length(xs) if !(xs[i] in elim)]

#   if length(xnew_str) == 0
#     phi = hom(R, cR, [cR(0) for i in 1:length(xs)])
#     ambR = codomain(phi)

#     if length(Igens) == 0
#       Inew = ideal(ambR, [ambR(0)])
#     else
#       Inew = ideal(ambR, phi.(Igens))
#     end
#     normal_Sgens = phi.(Sgens)
#     # This is a hack as phi.(Sgens) returns a list without types on the empty list
#     if normal_Sgens == Any[]
#       normal_Sgens = Vector{RingElem}()
#     end
#   else
#     Rnew, xnew = polynomial_ring(coefficient_ring(R), length(xnew_str))

#     zero_elim_var = elem_type(Rnew)[]
#     j = 1
#     for i in 1:length(xs)
#       if xs[i] in elim
#         push!(zero_elim_var, Rnew(0))
#       else
#         push!(zero_elim_var, xnew[j])
#         j += 1
#       end
#     end

#     phi = hom(R, Rnew, zero_elim_var)

#     ambR = codomain(phi)
#     if length(Igens) == 0
#       Inew = ideal(ambR, ambR(0))
#     else
#       Inew = ideal(ambR, phi.(Igens))
#     end

#     if length(Sgens) == 0
#       normal_Sgens = Vector{RingElem}()
#     else
#       Sgens_new = phi.(Sgens)
#       normal_Sgens = [normal_form(g, Inew) for g in Sgens_new]
#       if !(ambR(0) in normal_Sgens)
#         normal_Sgens = gens_2_prime_divisors(Sgens_new)
#       end
#     end
#   end

#   if isone(Inew) || ambR(0) in normal_Sgens
#     MRS_new = MatroidRealizationSpace(Inew, normal_Sgens, ambR, nothing, MRS.char, MRS.q, MRS.ground_ring)
#     set_attribute!(MRS_new, :is_realizable_sp, :false)
#     return MRS_new
#   end

#   Xnew = matrix(ambR, [phi(X[i, j]) for i in 1:nr, j in 1:nc])

#   #Try to reduce the matrix one last time using the ideal and the inequations_sp
#   m, n = size(Xnew)
#   A = base_ring(R)
#   if length(gens(Inew)) > 0 && gens(Inew)[1] != 0
#     for i in 1:m, j in 1:n
#       if A == ZZ
#         Xnew[i, j] = mod(Xnew[i, j], gens(Inew)[1])
#       else
#         Xnew[i, j] = reduce(Xnew[i, j], gens(Inew))
#       end
#     end
#   end

#   for j in 1:n
#     g = gcd(Xnew[:, j]...)
#     prime_divisors = poly_2_prime_divisors(g)
#     for f in prime_divisors
#       if f in normal_Sgens
#         for i in 1:m
#           Xnew[i, j] = Xnew[i, j] / f
#         end
#       end
#     end
#   end
#   MRS_new = MatroidRealizationSpace(Inew, normal_Sgens, ambR, Xnew, MRS.char, MRS.q, MRS.ground_ring)
#   return MRS_new
# end

# Exports
export MatroidRealizationSpace_SelfProj
export inequations
export is_realizable
export realization
export realization_matrix
export realization_space
export underlying_scheme
export selfproj_realization_ideal