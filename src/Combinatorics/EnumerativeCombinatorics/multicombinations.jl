@doc raw"""
    multicombinations(n::IntegerUnion, k::IntegerUnion)

Return an iterator over all $k$-combinations of ${1,...,n}$ with repetition,
produced in lexicographically ascending order.

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
multicombinations(n::T, k::T) where T<:IntegerUnion = MultiCombinations(Base.oneto(n), n, k)


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

MultiCombinations(v::AbstractArray, k::T) where {T<:IntegerUnion} = MultiCombinations(v, T(length(v)), k)


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
  return Combination{T}(C.v[state]), state
end

Base.length(C::MultiCombinations) = binomial(Int(C.n)+Int(C.k)-1, Int(C.k))

Base.eltype(::Type{<:MultiCombinations{T}}) where {T} = Combination{eltype(T)}

function Base.show(io::IO, C::MultiCombinations{<:Base.OneTo})
  print(io, "Iterator over the ", C.k, "-combinations of ", 1:C.n, " with repetition")
end

function Base.show(io::IO, C::MultiCombinations)
  print(io, "Iterator over the ", C.k, "-combinations of ", C.v, " with repetition")
end
