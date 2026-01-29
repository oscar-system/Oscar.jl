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
  betti_numbers = get_attribute!(
    () -> [ZZ(-1) for _ in 1:(d + 1)], v, :betti_number
  )::Vector{ZZRingElem}

  # compute the Betti number if needed
  k = i >> 1 # i is even, so divide by two and use that as index
  if betti_numbers[k + 1] == -1
    f_vector::Vector{Int} = pm_object(v).F_VECTOR
    pushfirst!(f_vector, 1)
    betti_numbers[k + 1] = ZZRingElem(
      sum((-1)^(i - k) * binomial(i, k) * f_vector[d - i + 1] for i in k:d)
    )
  end

  # return result
  return betti_numbers[k + 1]
end
