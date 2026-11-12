@doc raw"""
    det(
      c::AbsHyperComplex{T};
      lower_bound::Union{Int, Nothing},
      upper_bound::Union{Int, Nothing},
      direction::Symbol, reduction_map
    ) where {T<:OFPModule}

Compute the determinant of a one-dimensional complex as outlined for instance in [GKZ08](@cite). The keyword `direction` can be set to either `:from_left_to_right` to start from an upper bound on the non-zero modules, or `:from_right_to_left` to start from a lower bound. All modules in the complex need to be free and defined over an integral domain for which the field of fractions exists. 

The user may specify a `reduction_map` which is used to find suitable minors for the computation of the determinant. In the case of a polynomial 
ring this can, for instance, be the evaluation of the matrices at a random point and/or reduction modulo primes. 

!!! note
    It is assumed, but not checked that there exists a single range of integers `r = a:b` such that the modules `c[i]` are non-zero for `i` in `r` and zero otherwise. 
"""
function det(
    c::AbsHyperComplex{T};
    lower_bound::Union{Int, Nothing}=has_lower_bound(c, 1) ? lower_bound(c, 1) : nothing,
    upper_bound::Union{Int, Nothing}=has_upper_bound(c, 1) ? upper_bound(c, 1) : nothing,
    direction::Symbol= isnothing(lower_bound) ? :from_right_to_left : :from_left_to_right, 
    reduction_map=_generic_reduction_map(c, upper_bound, lower_bound)
  ) where {T <: OFPModule}
  @assert is_one(dim(c)) "complex must be 1-dimensional"
  # We have a `direction` of the complex, depending on whether this is a chain, or a 
  # cochain complex. In addition, we have a direction for the computations to be 
  # carried out, i.e. with which morphism to start and in which direction to proceed. 
  # While the first is predetermined by the input, there might be preferable choices 
  # for the latter. 
  return _det(c, Val{direction}, Val{Oscar.direction(c, 1)}, upper_bound, lower_bound; reduction_map)
end

# We try to provide reasonable automatic reduction maps. 
# The user can overwrite this for their use-cases, or pass their 
# own reduction map. False results can not be returned, as this 
# is only used to search for suitable minors, but their correctness
# is verified a posteriori.
function _generic_reduction_map(c::AbsHyperComplex, ::Nothing, ::Nothing)
  error("no upper or lower bound provided; can not find starting point to extract the ring")
end

function _generic_reduction_map(c::AbsHyperComplex, ub::Int, ::Any)
  return _generic_reduction_map(base_ring(c[ub]))
end

function _generic_reduction_map(c::AbsHyperComplex, ::Nothing, lb::Int)
  return _generic_reduction_map(base_ring(c[lb]))
end

# The default is no reduction map at all
_generic_reduction_map(::Ring) = nothing

# For polynomial rings over a field we can do better
function _generic_reduction_map(R::MPolyRing{<:FieldElem})
  kk = coefficient_ring(R)
  pt = elem_type(kk)[rand(kk, -100:100) for _ in 1:ngens(R)]
  return hom(R, kk, pt)
end

# For polynomials over the rationals we reduce modulo primes.
# This can (in theory!) run into divisions by zero. In that 
# case the user has to try again or use their own reduction map. 
# In practice, this is most likely never encountered. 
function _generic_reduction_map(R::MPolyRing{QQFieldElem})
  p = next_prime(2^50 + rand(0:2^10))
  kk = GF(p)
  pt = elem_type(kk)[rand(kk) for _ in 1:ngens(R)]
  return hom(R, kk, f->kk(numerator(f))*inv(kk(denominator(f))), pt)
end

# internal implementations depending on the constellation of 
# direction and type of the complex
#
# The first four methods are used for preparation of the arguments
# for the fifth method below. This is to make avoid `if`-`else`-
# constructions and make use of the dispatch. 
function _det(c::AbsHyperComplex, ::Type{Val{:from_left_to_right}}, 
    ::Type{Val{:cochain}}, 
    upper_bound::Union{Nothing, Int}, lower_bound::Int; 
    reduction_map=nothing
  )
  return _det(c, lower_bound, 1, 0; reduction_map)
end

function _det(c::AbsHyperComplex, ::Type{Val{:from_left_to_right}}, 
    ::Type{Val{:chain}},
    upper_bound::Union{Nothing, Int}, lower_bound::Int; 
    reduction_map=nothing
  )
  return _det(c, lower_bound, 1, 1; reduction_map)
end

function _det(c::AbsHyperComplex, ::Type{Val{:from_right_to_left}}, 
    ::Type{Val{:chain}}, 
    upper_bound::Int, lower_bound::Union{Nothing, Int};
    reduction_map=nothing
  )
  return _det(c, upper_bound, -1, 0; reduction_map)
end

function _det(c::AbsHyperComplex, ::Type{Val{:from_right_to_left}}, 
    ::Type{Val{:cochain}}, 
    upper_bound::Int, lower_bound::Union{Nothing, Int};
    reduction_map=nothing
  )
  return _det(c, upper_bound, -1, -1; reduction_map)
end

function _det(
    c::AbsHyperComplex, 
    ind::Int, # starting index and running variable
    increment::Int, # the direction in which to proceed
    offset::Int; # offset in the index for getting the relevant map
    reduction_map=nothing
  )
  while !can_compute_index(c, ind) || is_zero(c[ind]) 
    ind += increment
  end
  c[ind]::FreeMod
  R = base_ring(c[ind])
  result = one(fraction_field(R))
  r = ngens(c[ind]) # the rank of the current map
  I = collect(1:r)
  while can_compute_map(c, ind + offset) && !is_zero(map(c, ind + offset))
    A = !is_zero(offset) ? transpose(matrix(map(c, ind + offset))) : matrix(map(c, ind))
    J, p = _find_minor(A[I, :], reduction_map)
    # update variables
    if !is_zero(offset)
      result = is_even(ind) ? result*p : result//p
    else
      # outgoing maps
      result = is_even(ind) ? result//p : result*p
    end
    I = Int[i for i in 1:ncols(A) if !(i in J)]
    r = ncols(A) - r
    ind += increment
  end
  return result
end

# internal method to determine a non-zero maximal minor and return 
# a pair `(J, p)` of the multi-index `J` of that minor `p`.
function _find_minor(A::MatrixElem, reduction_map::Nothing)
  r = nrows(A)
  n = ncols(A)
  R = base_ring(A)
  for J in combinations(n, r)
    p = det(A[:, data(J)])
    is_zero(p) && continue
    return data(J), p
  end
  error("no non-zero minor of suitable size found; complex is not generically exact")
end

# the case of a non-trivial reduction map
function _find_minor(A::MatrixElem, reduction_map)
  r = nrows(A)
  n = ncols(A)
  R = base_ring(A)
  @assert reduction_map(one(R)) isa FieldElem "reduction map must produce elements in a field"
  A_red = map_entries(reduction_map, A)::MatrixElem{<:FieldElem}
  rk, M = rref(A_red)
  @assert r == rk "reduction map did not recover the anticipated rank ($r anticipated vs. $rk real)"
  pivots = Int[findfirst(!is_zero(a) for a in M[i, :]) for i in 1:nrows(M)]
  @assert length(pivots) == r "collection of pivot indices failed"
  p = det(A[:, pivots])
  @assert !is_zero(p) "the minor found via the reduction map was zero; try using another or no reduction map"
  return pivots, p
end

