########################################################################
# Ordered multiindices of the form 0 < i₁ < i₂ < … < iₚ ≤ n.
#
# For instance, these provide a coherent way to enumerate the generators 
# of an exterior power ⋀ ᵖ M of a module M given a chosen set of 
# generators of M. But they can also be used for other purposes 
# like enumerating subsets of p elements out of a set with n elements.
########################################################################
mutable struct OrderedMultiIndex{IntType<:IntegerUnion}
  i::Vector{IntType}
  n::IntType

  function OrderedMultiIndex(i::Vector{T}, n::T) where {T<:IntegerUnion}
    @assert all(k->i[k]<i[k+1], 1:length(i)-1) "indices must be strictly ordered"
    return new{T}(i, n)
  end
end

function ordered_multi_index(i::Vector{T}, n::T) where {T<:IntegerUnion} 
  return OrderedMultiIndex(i, n)
end

# For an ordered multiindex i = (0 < i₁ < i₂ < … < iₚ ≤ n) this returns the vector (i₁,…,iₚ).
indices(a::OrderedMultiIndex) = a.i

# For an ordered multiindex i = (0 < i₁ < i₂ < … < iₚ ≤ n) this returns n.
bound(a::OrderedMultiIndex) = a.n

# For an ordered multiindex i = (0 < i₁ < i₂ < … < iₚ ≤ n) this returns p.
length(a::OrderedMultiIndex) = length(a.i)

# For an ordered multiindex i = (0 < i₁ < i₂ < … < iₚ ≤ n) this returns iₖ.
getindex(a::OrderedMultiIndex, k::Int) = a.i[k]

index_type(a::OrderedMultiIndex) = index_type(typeof(a))
index_type(::Type{OrderedMultiIndex{T}}) where {T} = T

# Internal function for "multiplication" of ordered multiindices.
# 
# For i = (0 < i₁ < i₂ < … < iₚ ≤ n) and j = (0 < j₁ < j₂ < … < jᵣ ≤ n)
# the result is a pair `(sign, a)` with `sign` either 0 in case that 
# iₖ = jₗ for some k and l, or ±1 depending on the number of transpositions 
# needed to put (i₁, …, iₚ, j₁, …, jᵣ) into a strictly increasing order 
# to produce `a`.
function _mult(a::OrderedMultiIndex{T}, b::OrderedMultiIndex{T}) where {T}
  @assert bound(a) == bound(b) "multiindices must have the same bounds"

  # in case of a double index return zero
  any(in(indices(b)), indices(a)) && return 0, indices(a)

  p = length(a)
  q = length(b)
  result_indices = vcat(indices(a), indices(b))
  sign = 1

  # bubble sort result_indices and keep track of the sign
  for k in p:-1:1
    l = k
    c = result_indices[l]
    while l < p + q && c > result_indices[l+1]
      result_indices[l] = result_indices[l+1]
      sign = -sign
      l = l+1
    end
    result_indices[l] = c
  end
  return sign, result_indices
end

# For two ordered multiindices i = (0 < i₁ < i₂ < … < iₚ ≤ n) 
# and j = (0 < j₁ < j₂ < … < jᵣ ≤ n) this returns a pair `(sign, a)` 
# with `sign` either 0 in case that iₖ = jₗ for some k and l, 
# or ±1 depending on the number of transpositions needed to put 
# (i₁, …, iₚ, j₁, …, jᵣ) into a strictly increasing order to produce `a`.
function _wedge(a::OrderedMultiIndex{T}, b::OrderedMultiIndex{T}) where {T}
  sign, ind = _mult(a, b)
  iszero(sign) && return sign, a
  return sign, OrderedMultiIndex(ind, bound(a))
end

function _wedge(a::Vector{T}) where {T <: OrderedMultiIndex}
  isempty(a) && error("list must not be empty")
  isone(length(a)) && return 1, first(a)
  k = div(length(a), 2)
  b = a[1:k]
  c = a[k+1:end]
  sign_b, ind_b = _wedge(b)
  sign_c, ind_c = _wedge(c)
  sign, ind = _wedge(ind_b, ind_c)
  return sign * sign_b * sign_c, ind
end

function ==(a::OrderedMultiIndex{T}, b::OrderedMultiIndex{T}) where {T}
  return bound(a) == bound(b) && indices(a) == indices(b)
end

function Base.hash(a::OrderedMultiIndex, h::UInt)
  return hash(indices(a), hash(bound(a), h))
end

########################################################################
# A data type to facilitate iteration over all ordered multiindices 
# of the form 0 < i₁ < i₂ < … < iₚ ≤ n for fixed 0 ≤ p ≤ n.
#
# Example:
#
#   for i in OrderedMultiIndexSet(3, 5)
#      # do something with i = (0 < i₁ < i₂ < i₃ ≤ 5).
#   end
########################################################################
mutable struct OrderedMultiIndexSet
  n::Int
  p::Int

  function OrderedMultiIndexSet(p::Int, n::Int)
    @assert 0 <= p <= n "invalid bounds"
    return new(n, p)
  end
end

bound(I::OrderedMultiIndexSet) = I.n
index_length(I::OrderedMultiIndexSet) = I.p

Base.eltype(I::OrderedMultiIndexSet) = OrderedMultiIndex
Base.length(I::OrderedMultiIndexSet) = binomial(bound(I), index_length(I))

function Base.iterate(I::OrderedMultiIndexSet)
  ind = OrderedMultiIndex([i for i in 1:index_length(I)], bound(I))
  return ind, ind
end

function Base.iterate(I::OrderedMultiIndexSet, state::OrderedMultiIndex)
  bound(I) == bound(state) || error("index not compatible with set")
  ind = copy(indices(state))
  l = length(state)
  while l > 0 && ind[l] == bound(I) - length(state) + l
    l = l - 1
  end
  iszero(l) && return nothing
  ind[l] = ind[l] + 1
  l = l + 1
  while l <= length(state)
    ind[l] = ind[l-1] + 1
    l = l + 1
  end
  result = OrderedMultiIndex(ind, bound(I))
  return result, result
end

function Base.show(io::IO, ind::OrderedMultiIndex)
  i = indices(ind)
  print(io, "0 ")
  for i in indices(ind)
    print(io, "< $i ")
  end
  print(io, "<= $(bound(ind))")
end

# For an ordered multiindex i = (0 < i₁ < i₂ < … < iₚ ≤ n) this 
# returns the number k so that i appears at the k-th spot in the 
# enumeration of all ordered multiindices for this pair 0 ≤ p ≤ n.
function linear_index(ind::OrderedMultiIndex)
  n = bound(ind)
  p = length(ind)
  iszero(p) && return 1
  isone(p) && return ind[1]
  i = indices(ind)
  return binomial(n, p) - binomial(n - first(i) + 1, p) + linear_index(OrderedMultiIndex(i[2:end].-first(i), n-first(i)))
end

# For a pair 0 ≤ p ≤ n return the k-th ordered multiindex in the 
# enumeration of all ordered multiindices (0 < i₁ < i₂ < … < iₚ ≤ n).
function ordered_multi_index(k::Int, p::Int, n::Int)
  (k < 1 || k > binomial(n, p)) && error("index out of range")
  iszero(p) && return OrderedMultiIndex(Int[], n)
  isone(p) && return OrderedMultiIndex([k], n)
  n == p && return OrderedMultiIndex([k for k in 1:n], n)
  bin = binomial(n, p)
  i1 = findfirst(j->(bin - binomial(n - j, p) > k - 1), 1:n-p)
  if i1 === nothing 
    prev_res = ordered_multi_index(k - bin + 1, p-1, p-1)
    k = n-p+1
    return OrderedMultiIndex(pushfirst!(indices(prev_res).+k, k), n)
  else
    prev_res = ordered_multi_index(k - bin + binomial(n - i1 + 1, p), p-1, n - i1)
    return OrderedMultiIndex(pushfirst!(indices(prev_res).+i1, i1), n)
  end
end
