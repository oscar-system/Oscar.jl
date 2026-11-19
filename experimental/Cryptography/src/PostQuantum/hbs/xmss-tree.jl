# ---------------------------------------------------------------------------
# Merkle trees over Winternitz keys, the building block of both kinds of
# schemes. The tree number `tree` in the layer `layer` has the public keys of
# the Winternitz keys `(layer, tree, 0)`, ..., `(layer, tree, 2^hp - 1)` as leaves.
# ---------------------------------------------------------------------------

# The node hash of the tree number `tree` in the layer `layer`.
_xmss_node(n::Int, seed::Vector{UInt8}, layer::Int, tree::Int) =
    (height, index, left, right) -> _hbs_thash(n, seed, _hbs_adrs(_HBS_TREE, layer, tree, height, index), left, right)

# The tree of height `hp`; its construction needs the public keys of all Winternitz keys.
function _xmss_tree(n::Int, w::Int, hp::Int, sk_seed::Vector{UInt8}, seed::Vector{UInt8}, layer::Int, tree::Int)
    leaves = [_wots_pk(n, w, sk_seed, seed, (layer, tree, j)) for j in 0:2^hp - 1]
    return merkle_tree(leaves; node = _xmss_node(n, seed, layer, tree))
end

# Sign the message `m` of `n` bytes with the Winternitz key number `leaf` of the tree
# `levels`. Return the signature of `len * n` bytes followed by the authentication path.
function _xmss_sign(n::Int, w::Int, levels::Vector{Vector{Vector{UInt8}}}, m::Vector{UInt8}, sk_seed::Vector{UInt8}, seed::Vector{UInt8}, layer::Int, tree::Int, leaf::Int)
    sig = _wots_sign(n, w, m, sk_seed, seed, (layer, tree, leaf))
    return vcat(reduce(vcat, sig), reduce(vcat, merkle_auth_path(levels, leaf)))
end

# The root of the tree that a signature `sig` of `m` by the key number `leaf` belongs to.
function _xmss_root_from_sig(n::Int, w::Int, hp::Int, sig::Vector{UInt8}, m::Vector{UInt8}, seed::Vector{UInt8}, layer::Int, tree::Int, leaf::Int)
    len = _wots_len(n, w)
    @req length(sig) == (len + hp) * n "signature must have length $((len + hp) * n), got $(length(sig))"
    parts = _hbs_split(sig, n)
    pk = _wots_pk_from_sig(n, w, parts[1:len], m, seed, (layer, tree, leaf))
    return _merkle_root_from_path(pk, leaf, parts[len + 1:end]; node = _xmss_node(n, seed, layer, tree))
end
