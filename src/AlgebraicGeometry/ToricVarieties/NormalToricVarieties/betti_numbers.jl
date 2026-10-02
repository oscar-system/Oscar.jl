############################
# Betti numbers
############################

@doc raw"""
    betti_number(v::NormalToricVarietyType, i::Int)

Compute the `i`-th Betti number of the normal toric variety `v`. 
Specifically, this method returns the dimension of the i-th 
simplicial homology group (with rational coefficients) of `v`. 
The employed algorithm is derived from theorem 12.3.12 in 
[CLS11](@cite). Note that this theorem requires that the normal 
toric variety `v` is both complete and simplicial.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> betti_number(P3,0)
1

julia> betti_number(P3, 1)
0
```
"""
function betti_number(v::NormalToricVarietyType, i::Int)
  @req is_complete(v) && is_simplicial(v) "Currently, the computation of Betti numbers is limited to complete and simplicial toric varieties"

  # check input
  d = dim(v)::Int
  if !(0 <= i <= 2 * d) || isodd(i)
    return ZZRingElem(0)
  end

  # extract vector of currently-known Betti numbers (or create it if necessary)
  cached_betti_numbers = get_attribute!(
    () -> [ZZ(-1) for _ in 1:(d + 1)], v, :betti_numbers
  )::Vector{ZZRingElem}

  # compute the Betti number if needed
  k = i >> 1 # i is even, so divide by two and use that as index
  if cached_betti_numbers[k + 1] == -1
    f_vector::Vector{Int} = pm_object(v).F_VECTOR
    pushfirst!(f_vector, 1)
    cached_betti_numbers[k + 1] = ZZRingElem(
      sum((-1)^(i - k) * binomial(i, k) * f_vector[d - i + 1] for i in k:d)
    )
  end

  # return result
  return deepcopy(cached_betti_numbers[k + 1])
end

@doc raw"""
    betti_numbers(v::NormalToricVarietyType) -> Vector{ZZRingElem}

Return all ordinary rational Betti numbers of `v`, from ``b_0`` through
``b_{2d}``, where ``d`` is the dimension of `v`. The returned vector includes
the zero odd-degree entries.

The variety `v` must be complete and simplicial, but it need not be smooth. Use
[`betti_number`](@ref) to request a single degree.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2);

julia> betti_numbers(P2)
5-element Vector{ZZRingElem}:
 1
 0
 1
 0
 1
```
"""
function betti_numbers(v::NormalToricVarietyType)
  @req is_complete(v) && is_simplicial(v) "Currently, the computation of Betti numbers is limited to complete and simplicial toric varieties"
  return [betti_number(v, i) for i in 0:(2 * dim(v))]
end
