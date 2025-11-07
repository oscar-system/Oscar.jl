@doc raw"""
    combinations(n::IntegerUnion, k::IntegerUnion)

Return an iterator over all $k$-combinations of ${1,...,n}$, produced in
lexicographically ascending order.

# Examples

```jldoctest
julia> C = combinations(4, 3)
Iterator over the 3-combinations of 1:4

julia> collect(C)
4-element Vector{Combination{Int64}}:
 [1, 2, 3]
 [1, 2, 4]
 [1, 3, 4]
 [2, 3, 4]
```
"""
combinations(n::T, k::T; inplace::Bool = false) where T<:IntegerUnion = Combinations(Base.oneto(n), n, k, inplace)

@doc raw"""
    combinations(v::AbstractVector, k::IntegerUnion)

Return an iterator over all `k`-combinations of a given vector `v` produced in
lexicographically ascending order of the indices.

# Examples

```jldoctest
julia> C = combinations(['a', 'b', 'c', 'd'], 3)
Iterator over the 3-combinations of ['a', 'b', 'c', 'd']

julia> collect(C)
4-element Vector{Combination{Char}}:
 ['a', 'b', 'c']
 ['a', 'b', 'd']
 ['a', 'c', 'd']
 ['b', 'c', 'd']
```
"""
combinations(v::AbstractVector, k::IntegerUnion) = Combinations(v, k)

@inline function Base.iterate(C::Combinations{<:AbstractVector{T}, U}, state::Vector{U} = U[min(C.k - 1, i) for i in Base.oneto(C.k)]) where {T, U<:IntegerUnion}
  if is_zero(C.k) # special case to generate 1 result for k = 0
    if isempty(state)
      return Combination{T}(T[]), U[0]
    end
    return nothing
  end
  for i in C.k:U(-1):U(1)
    state[i] += 1
    if state[i] > (C.n - (C.k - i))
      continue
    end
    for j in i+1:C.k
      state[j] = state[j - 1] + 1
    end
    break
  end
  if state[1] > C.n - C.k + 1
    return nothing
  end
  return Combination{T}(C.v[state]), state
end

Base.length(C::Combinations) = binomial(Int(C.n), Int(C.k))

Base.eltype(::Type{<:Combinations{T}}) where {T} = Combination{eltype(T)}

function Base.show(io::IO, C::Combinations{<:Base.OneTo})
  print(io, "Iterator over the ", C.k, "-combinations of ", 1:C.n)
end

function Base.show(io::IO, C::Combinations)
  print(io, "Iterator over the ", C.k, "-combinations of ", C.v)
end

################################################################################
#
#  Array-like functionality for Combination
#
################################################################################

function Base.show(io::IO, ::MIME"text/plain", P::Combination)
  p = data(P)
  if isempty(p)
    print(io, "Empty combination")
    return
  end
  print(io, p)
end

data(C::Combination) = C.v

function Base.size(C::Combination)
  return size(data(C))
end

function Base.length(C::Combination)
  return length(data(C))
end

function Base.getindex(C::Combination, i::IntegerUnion)
  return getindex(data(C), Int(i))
end

function Base.copy(C::Combination{T}) where T
  return Combination{T}(copy(data(C)))
end

function Base.getindex(C::Combinations{Base.OneTo}, i::IntegerUnion)
  return combination(C.n, C.k, i)
end

function Base.getindex(C::Combinations, i::IntegerUnion)
  c = combination(C.n, C.k, i)
  return C.v[data(c)]
end

################################################################################
#
#  Ranking / unranking combinations
#
################################################################################


@doc raw"""
    linear_index(C::Combination, n::IntegerUnion)

For a combination `C` containing `k` elements from `1:n`,
return the index `i` so that `C` occurs as the `i`th element in the
iteration over all such combinations (in lexicographic order).

The algorithm is as described in [Lie97; Section 3.3](@cite).

# Examples

```jldoctest
julia> C = Combination([1, 3, 5, 7])
[1, 3, 5, 7]

julia> Oscar.linear_index(C, 7)
15
```
"""
function linear_index(C::Combination, n::IntegerUnion)
  @req is_empty(C) || C[end] <= n "Combination must be bounded by n"
  k = length(C)
  iszero(k) && return 1
  isone(k) && return Int(C[1])

  r = binomial(n,k)
  i = 1
  while i <= k
    r -= binomial(n-C[i], k-i+1)
    i += 1
    k-i+1 == n-C[i-1] && return r
  end
  return r
end

@doc raw"""
    combination(n::IntegerUnion, k::IntegerUnion, r::IntegerUnion)

Return the `r`th combination in the iteration over `combinations(n,k)`.

The algorithm is as described in [Lie97; Section 3.3](@cite).

# Examples

```jldoctest
julia> C1 = collect(combinations(15, 3))[13]
[1, 2, 15]

julia> C1 = Oscar.combination(15, 3, 13)
[1, 2, 15]
```
"""
function combination(n::T, k::IntegerUnion, r::IntegerUnion) where T<:IntegerUnion
  @req 1 <= r <= binomial(Int(n), Int(k)) "index out of range"
  iszero(k) && return Combination{T}(T[])
  isone(k) && return Combination{T}(T[r])
  n == k && return Combination{T}(collect(T(1):n))

  C = zeros(T, k)
  r = binomial(Int(n),Int(k)) - r
  j = T(1)
  for i in 1:k
    b = binomial(Int(n) - j, Int(k) - i + 1)
    while b > r
      j += T(1)
      b = binomial(Int(n) - j, Int(k) - i + 1)
    end
    C[i] = j
    j += T(1)
    r -= b
  end
  return Combination{T}(C)
end


################################################################################
#
#  Misc (for OrderedMultiIndex compat)
#
################################################################################

# merge sort vcat(a,b) and keep track of the sign of permutation.
# requires that a and b are both sorted, and contain no common elements.
function merge_sorted_with_sign(a::Vector{T}, b::Vector{T}) where T<:IntegerUnion
  isdisjoint(a,b) || return nothing, 0

  result = zeros(T, length(a)+length(b))
  p = length(a)
  q = length(b)

  i = 1
  j = 1
  sign = 1
  while i <= p || j <= q
    if i > p
      result[i+j-1] = b[j]
      j += 1
    elseif j > q || a[i] < b[j]
      result[i+j-1] = a[i]
      i += 1
    else
      result[i+j-1] = b[j]
      is_even(p-i) && (sign = -sign)
      j += 1
    end
  end

  return result, sign
end


# For two combinations a = [i₁, i₂, …, iₚ] and b = [j₁, j₂, …, jᵣ]
# the result is a pair `(sign, c)` with `c` a new combination,
# and `sign` either 0 in the case that
# iₖ = jₗ for some k and l, or ±1 depending on the number of transpositions
# needed to put [i₁, …, iₚ, j₁, …, jᵣ] into a strictly increasing order
# to produce `c`. If the combinations are not disjoint,
# then we return `nothing` for `c`.
function _wedge(a::Combination{T}, b::Combination{T}) where {T}
  c, sign = merge_sorted_with_sign(data(a), data(b))
  sign == 0 && return sign, c
  return sign, Combination{T}(c)
end


# # This could be optimized, but don't think it is currently used anywhere so leaving it for now.
function _wedge(a::Vector{T}) where {T <: Combination}
  @req !isempty(a) "list must not be empty"
  isone(length(a)) && return 1, first(a)
  k = div(length(a), 2)
  b = a[1:k]
  c = a[k+1:end]
  sign_b, comb_b = _wedge(b)
  sign_c, comb_c = _wedge(c)
  sign, comb = _wedge(comb_b, comb_c)
  return sign * sign_b * sign_c, comb
end
