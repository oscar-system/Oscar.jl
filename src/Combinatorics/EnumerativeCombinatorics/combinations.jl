@doc raw"""
    combinations(n::Int, k::Int)

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
combinations(n::Int, k::Int) = Combinations(1:n, k)

@doc raw"""
    combinations(v::AbstractVector, k::Int)

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
function combinations(v::AbstractVector{T}, k::Int) where T
  return Combinations(v, k)
end

Combinations(v::AbstractArray{T}, k::Int) where T = Combinations(v, length(v), k)

@inline function Base.iterate(C::Combinations, state = [min(C.k - 1, i) for i in 1:C.k])
  if C.k == 0 # special case to generate 1 result for k = 0
    if isempty(state)
      return Combination(eltype(C.v)[]), [0]
    end
    return nothing
  end
  for i in C.k:-1:1
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
  return Combination(C.v[state]), state
end

Base.length(C::Combinations) = binomial(C.n, C.k)

Base.eltype(::Type{Combinations{T}}) where {T} = Combination{eltype(T)}

function Base.show(io::IO, C::Combinations)
  print(io, "Iterator over the ", C.k, "-combinations of ", C.v)
end

################################################################################
#
#  Array-like functionality
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

data(P::Combination) = P.v

function Base.size(P::Combination)
  return size(data(P))
end

function Base.length(P::Combination)
  return length(data(P))
end

function Base.getindex(P::Combination, i::IntegerUnion)
  return getindex(data(P), Int(i))
end

function Base.setindex!(P::Combination, x::IntegerUnion, i::IntegerUnion)
  return setindex!(data(P), x, i)
end

function Base.copy(P::Combination)
  return Combination(copy(data(P)))
end
