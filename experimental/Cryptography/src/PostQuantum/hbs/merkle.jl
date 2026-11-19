# ---------------------------------------------------------------------------
# Merkle trees. A tree is stored as the vector of its levels: level 0 holds the
# leaves, the last level holds the root. The hash of two children is passed as
# the function `node(height, index, left, right)`, so that the schemes below
# can bind it to a position in the tree. By default it is plain SHAKE256.
# ---------------------------------------------------------------------------

_merkle_node(height::Int, index::Int, left::Vector{UInt8}, right::Vector{UInt8}) =
    shake_256(vcat(left, right), length(left))

@doc raw"""
    merkle_tree(leaves::Vector{Vector{UInt8}}; node = ...) -> Vector{Vector{Vector{UInt8}}}

Build the Merkle tree over `leaves`, whose number must be a power of two.
The result is the vector of the levels of the tree, starting with the leaves
and ending with the level that contains only the root. A node is the hash of
its two children, computed by `node(height, index, left, right)`; by default
this is the first `n` bytes of SHAKE256 applied to the concatenation, where
`n` is the length of `left`.

# Examples
```jldoctest
julia> leaves = [UInt8[i] for i in 1:4];

julia> levels = merkle_tree(leaves);

julia> length.(levels)
3-element Vector{Int64}:
 4
 2
 1

julia> merkle_root(levels) == merkle_root(merkle_tree(leaves))
true
```
"""
function merkle_tree(leaves::Vector{Vector{UInt8}}; node = _merkle_node)
    @req !isempty(leaves) && ispow2(length(leaves)) "the number of leaves must be a power of two, got $(length(leaves))"
    levels = [leaves]
    while length(levels[end]) > 1
        prev = levels[end]
        height = length(levels)
        push!(levels, [node(height, i - 1, prev[2i - 1], prev[2i]) for i in 1:length(prev) ÷ 2])
    end
    return levels
end

@doc raw"""
    merkle_root(levels::Vector{Vector{Vector{UInt8}}}) -> Vector{UInt8}

Return the root of the Merkle tree `levels` created by [`merkle_tree`](@ref).
"""
merkle_root(levels::Vector{Vector{Vector{UInt8}}}) = levels[end][1]

@doc raw"""
    merkle_auth_path(levels::Vector{Vector{Vector{UInt8}}}, index::Int) -> Vector{Vector{UInt8}}

Return the authentication path of the leaf with the zero-based position `index`
in the Merkle tree `levels`. It consists of the siblings of the nodes on the way
from the leaf to the root, starting at the leaf.

# Examples
```jldoctest
julia> levels = merkle_tree([UInt8[i] for i in 1:8]);

julia> length(merkle_auth_path(levels, 5))
3
```
"""
function merkle_auth_path(levels::Vector{Vector{Vector{UInt8}}}, index::Int)
    @req 0 <= index < length(levels[1]) "leaf index $index out of range"
    return [levels[h][(index >> (h - 1)) ⊻ 1 + 1] for h in 1:length(levels) - 1]
end

# Compute the root from the leaf at position `index` and its authentication path.
function _merkle_root_from_path(leaf::Vector{UInt8}, index::Int, path::Vector{Vector{UInt8}}; node = _merkle_node)
    x = leaf
    for (h, sibling) in enumerate(path)
        i = index >> (h - 1)
        x = iseven(i) ? node(h, i >> 1, x, sibling) : node(h, i >> 1, sibling, x)
    end
    return x
end

@doc raw"""
    merkle_verify(root::Vector{UInt8}, leaf::Vector{UInt8}, index::Int, path::Vector{Vector{UInt8}}; node = ...) -> Bool

Return `true` if `path` is the authentication path of `leaf` at the zero-based
position `index` in the Merkle tree with the given `root`. The function `node`
must be the one used to build the tree.

# Examples
```jldoctest
julia> levels = merkle_tree([UInt8[i] for i in 1:8]);

julia> path = merkle_auth_path(levels, 5);

julia> merkle_verify(merkle_root(levels), UInt8[6], 5, path)
true

julia> merkle_verify(merkle_root(levels), UInt8[6], 4, path)
false
```
"""
function merkle_verify(root::Vector{UInt8}, leaf::Vector{UInt8}, index::Int, path::Vector{Vector{UInt8}}; node = _merkle_node)
    return _merkle_root_from_path(leaf, index, path; node) == root
end
