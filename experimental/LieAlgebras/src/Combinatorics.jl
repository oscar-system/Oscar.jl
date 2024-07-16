# This file contains replacements for some functions in Combinatorics.jl.
# Most implementations here are quite slow and should be replaced by a
# more efficient implementation.

function combinations(n::Integer, k::Integer)
  return sort(subsets(n, k))
end

@doc raw"""
    combinations(v::AbstractVector{T}, k::Integer) where {T}

Return an iterator over all combinations of `k` elements of `v`.
In each iteration, the elements are returned in the order they appear in `v`.
The order of the combinations is lexicographic.

```jldoctest; setup = :(using Oscar.LieAlgebras)
julia> collect(combinations([1, 2, 3, 4], 2))
6-element Vector{Vector{Int64}}:
 [1, 2]
 [1, 3]
 [1, 4]
 [2, 3]
 [2, 4]
 [3, 4]
```
"""
function combinations(v::AbstractVector{T}, k::Integer) where {T}
  reorder(v, inds) = [v[i] for i in inds]
  return (reorder(v, inds) for inds in combinations(length(v), k))
end

function multicombinations(n::Integer, k::Integer)
  return sort(
    map(
      reverse ∘ Vector,
      reduce(vcat, collect(partitions(i, k, 1, n)) for i in 0:(k * n); init=Vector{Int}[]),
    ),
  )
end

@doc raw"""
    multicombinations(v::AbstractVector{T}, k::Integer) where {T}

Return an iterator over all combinations of `k` elements of `v` with repetitions.
In each iteration, the elements are returned in the order they appear in `v`.
The order of the combinations is lexicographic.

```jldoctest; setup = :(using Oscar.LieAlgebras)
julia> collect(multicombinations([1, 2, 3, 4], 2))
10-element Vector{Vector{Int64}}:
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
function multicombinations(v::AbstractVector{T}, k::Integer) where {T}
  reorder(v, inds) = [v[i] for i in inds]
  return (reorder(v, inds) for inds in multicombinations(length(v), k))
end

function permutations(n::Integer)
  return map(p -> Vector{Int}(p), AbstractAlgebra.SymmetricGroup(n))
end

@doc raw"""
    permutations(v::AbstractVector{T}) where {T}

Return an iterator over all permutations of `v`.
There is no guarantee on the order of the permutations.

```jldoctest; setup = :(using Oscar.LieAlgebras)
julia> sort(collect(permutations([1, 2, 3])))
6-element Vector{Vector{Int64}}:
 [1, 2, 3]
 [1, 3, 2]
 [2, 1, 3]
 [2, 3, 1]
 [3, 1, 2]
 [3, 2, 1]
```
"""
function permutations(v::AbstractVector{T}) where {T}
  reorder(v, p) = [v[p[i]] for i in 1:length(v)]
  return (reorder(v, p) for p in AbstractAlgebra.SymmetricGroup(length(v)))
end

function permutations_with_sign(n::Integer)
  return map(p -> (Vector{Int}(p), sign(p)), AbstractAlgebra.SymmetricGroup(n))
end

@doc raw"""
    permutations_with_sign(v::AbstractVector{T}) where {T}

Return an iterator over all permutations of `v` with their sign.
There is no guarantee on the order of the permutations.

```jldoctest; setup = :(using Oscar.LieAlgebras)
julia> sort(collect(permutations_with_sign([1, 2, 3])))
6-element Vector{Tuple{Vector{Int64}, Int64}}:
 ([1, 2, 3], 1)
 ([1, 3, 2], -1)
 ([2, 1, 3], -1)
 ([2, 3, 1], 1)
 ([3, 1, 2], 1)
 ([3, 2, 1], -1)
```
"""
function permutations_with_sign(v::AbstractVector{T}) where {T}
  reorder(v, p) = [v[p[i]] for i in 1:length(v)]
  return ((reorder(v, p), sign(p)) for p in AbstractAlgebra.SymmetricGroup(length(v)))
end
