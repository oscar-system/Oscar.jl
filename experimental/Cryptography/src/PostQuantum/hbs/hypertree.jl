# ---------------------------------------------------------------------------
# The stateless scheme with a hypertree of several Merkle trees, as in SPHINCS
# [BHH+15](@cite) and SPHINCS+ [BHK+19](@cite). There are `d` layers of Merkle
# trees of height `hp = h / d` over Winternitz keys. The top layer has one tree,
# whose root is part of the public key; each Winternitz key of a layer signs the
# root of a tree in the layer above. The Winternitz keys of the bottom layer sign
# the public keys of $2^h$ few-time keys (HORST or FORS). The digest of a message
# selects the few-time key, which signs a second part of the digest.
#
# Layout of the keys:
# - public key: `pk_seed || root`;
# - secret key: `sk_seed || sk_prf || pk_seed || root`;
# - signature: `R || few-time signature || d times (wots signature || path)`,
#   where `R` is the randomizer of the message.
#
# The few-time keys have the locations `(0, tree, leaf)`, the Winternitz keys
# in the layer `l` the locations `(l, tree, leaf)`.
# ---------------------------------------------------------------------------

# Number of bytes of the digest: the index of the few-time key and the part
# that is signed by it.
function _ht_md_sizes(p::HBS)
    fts_md = p.fts == :fors ? _fors_md_size(p.tau, p.k) : _horst_md_size(p.tau, p.k)
    return cld(p.h, 8), fts_md
end

# The few-time signature of the digest `md` with the key at `loc`.
function _ht_fts_sign(p::HBS, md::Vector{UInt8}, sk_seed::Vector{UInt8}, seed::Vector{UInt8}, loc)
    p.fts == :fors && return _fors_sign(p.n, p.tau, p.k, md, sk_seed, seed, loc)
    return _horst_sign(p.n, p.tau, p.k, p.x, md, sk_seed, seed, loc)
end

# The public key that a few-time signature of `md` belongs to (`nothing` for an invalid HORST signature).
function _ht_fts_pk_from_sig(p::HBS, sig::Vector{UInt8}, md::Vector{UInt8}, seed::Vector{UInt8}, loc)
    p.fts == :fors && return _fors_pk_from_sig(p.n, p.tau, p.k, sig, md, seed, loc)
    return _horst_pk_from_sig(p.n, p.tau, p.k, p.x, sig, md, seed, loc)
end

# Split the digest into the index of the few-time key (bytes in front) and the
# digest for the few-time signature, and the index into the tree and the leaf.
function _ht_split_digest(p::HBS, md::Vector{UInt8})
    isz, _ = _ht_md_sizes(p)
    idx = _hbs_indices(md[1:isz], 1, p.h)[1]
    return md[isz + 1:end], idx >> p.hp, idx & (2^p.hp - 1)
end

function _ht_keygen(p::HBS, rng::AbstractRNG)
    n = p.n
    sk_seed, sk_prf, pk_seed = rand(rng, UInt8, n), rand(rng, UInt8, n), rand(rng, UInt8, n)
    root = merkle_root(_xmss_tree(n, p.w, p.hp, sk_seed, pk_seed, p.d - 1, 0))
    return vcat(pk_seed, root), vcat(sk_seed, sk_prf, pk_seed, root)
end

function _ht_sign(p::HBS, msg::Vector{UInt8}, sk::Vector{UInt8})
    n = p.n
    @req length(sk) == p.sk_sz "secret key must have length $(p.sk_sz), got $(length(sk))"
    sk_seed, sk_prf, pk_seed, root = sk[1:n], sk[n + 1:2n], sk[2n + 1:3n], sk[3n + 1:4n]
    r = _hbs_prf_msg(n, sk_prf, msg)
    isz, fsz = _ht_md_sizes(p)
    fmd, tree, leaf = _ht_split_digest(p, _hbs_hmsg(isz + fsz, r, pk_seed, root, msg))
    sig = vcat(r, _ht_fts_sign(p, fmd, sk_seed, pk_seed, (0, tree, leaf)))
    m = _ht_fts_pk_from_sig(p, sig[n + 1:end], fmd, pk_seed, (0, tree, leaf))
    for layer in 0:p.d - 1
        levels = _xmss_tree(n, p.w, p.hp, sk_seed, pk_seed, layer, tree)
        append!(sig, _xmss_sign(n, p.w, levels, m, sk_seed, pk_seed, layer, tree, leaf))
        m = merkle_root(levels)
        leaf = tree & (2^p.hp - 1)
        tree >>= p.hp
    end
    return sig
end

function _ht_verify(p::HBS, sig::Vector{UInt8}, msg::Vector{UInt8}, pk::Vector{UInt8})
    n = p.n
    @req length(sig) == p.sig_sz "signature must have length $(p.sig_sz), got $(length(sig))"
    @req length(pk) == p.pk_sz "public key must have length $(p.pk_sz), got $(length(pk))"
    pk_seed, root = pk[1:n], pk[n + 1:2n]
    r = sig[1:n]
    isz, fsz = _ht_md_sizes(p)
    fmd, tree, leaf = _ht_split_digest(p, _hbs_hmsg(isz + fsz, r, pk_seed, root, msg))
    fts_sz = _fts_sig_size(p.fts, n, p.tau, p.k, p.x)
    m = _ht_fts_pk_from_sig(p, sig[n + 1:n + fts_sz], fmd, pk_seed, (0, tree, leaf))
    m === nothing && return false
    lsz = (p.len + p.hp) * n
    for layer in 0:p.d - 1
        part = sig[n + fts_sz + layer * lsz + 1:n + fts_sz + (layer + 1) * lsz]
        m = _xmss_root_from_sig(n, p.w, p.hp, part, m, pk_seed, layer, tree, leaf)
        leaf = tree & (2^p.hp - 1)
        tree >>= p.hp
    end
    return m == root
end
