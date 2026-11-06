# ---------------------------------------------------------------------------
# HORST (hash to obtain random subset with tree), the few-time signature scheme
# of SPHINCS [BHH+15](@cite). A key is one Merkle tree of height `tau` whose
# leaves are hashes of secret values. A digest is split into `k` indices of
# `tau` bits, and the signature reveals the secret values at these indices with
# their authentication paths. To make the paths shorter, the signature also
# contains the `2^x` nodes of the tree at height `tau - x`; the paths only lead
# up to these nodes and the root is computed from them. The public key is the root.
#
# A key is identified by `loc = (layer, tree, keypair)` as for Winternitz keys.
# ---------------------------------------------------------------------------

# Number of bytes of the digest that is signed.
_horst_md_size(tau::Int, k::Int) = cld(k * tau, 8)

# The tree of the key at `loc`. The secret values and the leaves are the
# same as for a FORS key with one tree, and so is the node hash.
function _horst_tree(n::Int, tau::Int, sk_seed::Vector{UInt8}, seed::Vector{UInt8}, loc)
    return _fors_tree(n, tau, sk_seed, seed, loc, 0)
end

# The public key of the key at `loc`.
_horst_pk(n::Int, tau::Int, sk_seed::Vector{UInt8}, seed::Vector{UInt8}, loc) =
    merkle_root(_horst_tree(n, tau, sk_seed, seed, loc))

# Sign the digest `md`; the signature has `k * (tau - x + 1) * n + 2^x * n` bytes.
function _horst_sign(n::Int, tau::Int, k::Int, x::Int, md::Vector{UInt8}, sk_seed::Vector{UInt8}, seed::Vector{UInt8}, loc)
    idx = _hbs_indices(md, k, tau)
    levels = _horst_tree(n, tau, sk_seed, seed, loc)
    sig = UInt8[]
    for i in idx
        append!(sig, _fors_sk(n, sk_seed, loc, i))
        append!(sig, reduce(vcat, merkle_auth_path(levels, i)[1:tau - x]; init = UInt8[]))
    end
    append!(sig, reduce(vcat, levels[tau - x + 1]))
    return sig
end

# The public key that a signature of the digest `md` belongs to, or `nothing`
# if a revealed value does not fit to the nodes at height `tau - x` in the signature.
function _horst_pk_from_sig(n::Int, tau::Int, k::Int, x::Int, sig::Vector{UInt8}, md::Vector{UInt8}, seed::Vector{UInt8}, loc)
    sz = (tau - x + 1) * n
    @req length(sig) == k * sz + 2^x * n "signature must have length $(k * sz + 2^x * n), got $(length(sig))"
    idx = _hbs_indices(md, k, tau)
    node = _fors_node(n, seed, loc, 0)
    top = _hbs_split(sig[k * sz + 1:end], n)
    for (j, i) in enumerate(idx)
        parts = _hbs_split(sig[(j - 1) * sz + 1:j * sz], n)
        leaf = _fors_leaf(n, seed, loc, i, parts[1])
        _merkle_root_from_path(leaf, i, parts[2:end]; node) == top[(i >> (tau - x)) + 1] || return nothing
    end
    return merkle_root(merkle_tree(top; node = (h, i, l, r) -> node(h + tau - x, i, l, r)))
end
