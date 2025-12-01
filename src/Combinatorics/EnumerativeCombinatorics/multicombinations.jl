@doc raw"""
    multicombinations(n::IntegerUnion, k::IntegerUnion; inplace::Bool=false)

Return an iterator over all $k$-combinations of ${1,...,n}$ with repetition,
produced in lexicographically ascending order.

If `inplace` is `true`, the elements of the iterator may share their memory. This
means that an element returned by the iterator may be overwritten 'in place' in
the next iteration step. This may result in significantly fewer memory allocations.
However, using the in-place version is only meaningful
if just one element of the iterator is needed at any time.
For example, calling `collect` on this iterator will not give useful results.

# Examples

```jldoctest
julia> C = multicombinations(4, 2)
Iterator over the 2-combinations of 1:4 with repetition

julia> collect(C)
10-element Vector{Combination{Int64}}:
 [1, 1]
 [1, 2]
 [1, 3]
 [1, 4]
 [2, 2]
 [2, 3]
 [2, 4]
 [3, 3]
 [3, 4]
 [4, 4]
```
"""
multicombinations(n::T, k::T; inplace::Bool = false) where T<:IntegerUnion = MultiCombinations(Base.oneto(n), n, k, inplace)


@doc raw"""
    multicombinations(v::AbstractVector{T}, k::Integer) where {T}

Return an iterator over all combinations of `k` elements of `v` with repetitions.
In each iteration, the elements are returned in the order they appear in `v`.
The order of the combinations is lexicographic.

```jldoctest
julia> collect(multicombinations(['a', 'b', 'c'], 2))
6-element Vector{Combination{Char}}:
 ['a', 'a']
 ['a', 'b']
 ['a', 'c']
 ['b', 'b']
 ['b', 'c']
 ['c', 'c']
```
"""
multicombinations(v::AbstractVector, k::IntegerUnion) = MultiCombinations(v, k)


@inline function Base.iterate(C::MultiCombinations{<:AbstractVector{T}, U}, state::Vector{U} = C.k == 0 ? U[] : (s = ones(U, C.k); s[Int(C.k)] = U(0); s)) where {T, U<:IntegerUnion}
  if is_zero(C.k) # special case to generate 1 result for k = 0
    if isempty(state)
      return Combination{T}(T[]), U[0]
    end
    return nothing
  end
  for i in C.k:U(-1):U(1)
    state[i] += 1
    if state[i] > C.n
      continue
    end
    for j in i+1:C.k
      state[j] = state[i]
    end
    break
  end
  if state[1] > C.n
    return nothing
  end

  c = C.inplace ? Combination{T}(state) : Combination{T}(C.v[state])
  return c, state
end

Base.eltype(::Type{<:MultiCombinations{T}}) where {T} = Combination{eltype(T)}

@doc raw"""
    number_of_multicombinations(n::IntegerUnion, k::IntegerUnion)

Return the number of $k$-combinations of ${1, \ldots, n}$.
If `n < 0` or `k < 0`, return `0`.
"""
function number_of_multicombinations(n::IntegerUnion, k::IntegerUnion)
  if n < 0 || k < 0
    return ZZ(0)
  end
  return binomial(ZZ(n) + ZZ(k) - 1, ZZ(k))
end

Base.length(C::MultiCombinations) = BigInt(number_of_multicombinations(C.n, C.k))

function Base.show(io::IO, C::MultiCombinations{<:Base.OneTo})
  print(io, "Iterator over the ", C.k, "-combinations of ", 1:C.n, " with repetition")
end

function Base.show(io::IO, C::MultiCombinations)
  print(io, "Iterator over the ", C.k, "-combinations of ", C.v, " with repetition")
end
