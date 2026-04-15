# This file contains replacements for some functions in Combinatorics.jl.
# Most implementations here are quite slow and should be replaced by a
# more efficient implementation.

function permutations(n::Integer)
  return map(Vector{Int}, AbstractAlgebra.SymmetricGroup(n))
end

@doc raw"""
    Oscar.LieAlgebras.permutations(v::AbstractVector{T}) where {T}

Return an iterator over all permutations of `v`.
There is no guarantee on the order of the permutations.

```jldoctest
julia> sort(collect(Oscar.LieAlgebras.permutations([1, 2, 3])))
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
    Oscar.LieAlgebras.permutations_with_sign(v::AbstractVector{T}) where {T}

Return an iterator over all permutations of `v` with their sign.
There is no guarantee on the order of the permutations.

```jldoctest
julia> sort(collect(Oscar.LieAlgebras.permutations_with_sign([1, 2, 3])))
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
