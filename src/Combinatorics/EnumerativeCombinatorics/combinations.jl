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

function Base.setindex!(C::Combination, x::IntegerUnion, i::IntegerUnion)
  return setindex!(data(C), x, i)
end

function Base.copy(C::Combination)
  return Combination(copy(data(C)))
end


################################################################################
#
#  Ranking / unranking combinations
#
################################################################################

# For a combination C containing k elements from 1:n,
# return the index i so that C occurs as the ith element in the
# iteration over all such combinations (which are lexicographically ordered).
# NOTE We use the method described in J. Liebehenschel - Ranking and unranking of lexicographically ordered words: an average-case analysis. Journal of Automata, Languages, and Combinatorics (1997).
# NOTE we have added an extra argument, since this info is not stored with
# a Combination as it is with an OrderedMultiIndex
function linear_index(C::Combination, n::IntegerUnion)
  k = length(C)
  iszero(k) && return 1
  isone(k) && return Int(C[1])

  r = binomial(n,k)
  i = 1
  while i < k+1
    r -= binomial(n-C[i], k-i+1)
    i += 1
    k-i+1 == n-C[i-1] && return r
  end
  return r
end

# Return the rth combination in the iteration over combinations(n,k)
# NOTE We use the method described in J. Liebehenschel - Ranking and unranking of lexicographically ordered words: an average-case analysis. Journal of Automata, Languages, and Combinatorics (1997).
# NOTE arguments are reordered to match combinatorial conventions
function combination(n::Int, k::Int, r::Int)
  (r < 1 || r > binomial(n, k)) && error("index out of range")
  iszero(k) && return Combination(Int[])
  isone(k) && return Combination([r])
  n == k && return Combination(collect(1:n))

  C = zeros(Int, k)
  r = binomial(n,k) - r
  j = 1
  for i in 1:k
    b = binomial(n - j, k - i + 1)
    while b > r - 1
      j += 1
      b = binomial(n - j, k - i _ 1)
    end
    C[i] = j
    j += 1
    r -= b
  end
  return Combination(C)
end
