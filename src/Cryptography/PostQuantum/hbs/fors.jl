# ---------------------------------------------------------------------------
# FORS (forest of random subsets), the few-time signature scheme of SPHINCS+
# [BHK+19](@cite). A key consists of `k` Merkle trees of height `a`, whose leaves
# are hashes of secret values. A digest is split into `k` indices of `a` bits,
# and the signature reveals the secret value at the index in each tree together
# with its authentication path. The public key is the compression of the roots.
# Signing more than a few messages with one key leaks more and more secret values.
#
# A key is identified by `loc = (layer, tree, keypair)` as for Winternitz keys.
# ---------------------------------------------------------------------------

# Number of bytes of the digest that is signed.
_fors_md_size(a::Int, k::Int) = cld(k * a, 8)

_fors_sk(n::Int, sk_seed::Vector{UInt8}, loc, leaf::Int) =
    _hbs_prf(n, sk_seed, _hbs_adrs(_HBS_PRF, _HBS_FTS_HASH, loc..., leaf))

_fors_leaf(n::Int, seed::Vector{UInt8}, loc, leaf::Int, sk::Vector{UInt8}) =
    _hbs_thash(n, seed, _hbs_adrs(_HBS_FTS_HASH, loc..., leaf), sk)

# The node hash inside of the tree number `i` of the key at `loc`.
_fors_node(n::Int, seed::Vector{UInt8}, loc, i::Int) =
    (height, index, left, right) -> _hbs_thash(n, seed, _hbs_adrs(_HBS_FTS_TREE, loc..., i, height, index), left, right)

# The tree number `i` (starting at 0) of the key at `loc`.
function _fors_tree(n::Int, a::Int, sk_seed::Vector{UInt8}, seed::Vector{UInt8}, loc, i::Int)
    leaves = [_fors_leaf(n, seed, loc, i * 2^a + j, _fors_sk(n, sk_seed, loc, i * 2^a + j)) for j in 0:2^a - 1]
    return merkle_tree(leaves; node = _fors_node(n, seed, loc, i))
end

# Compress the roots of the trees to the public key.
_fors_compress(n::Int, seed::Vector{UInt8}, loc, roots::Vector{Vector{UInt8}}) =
    _hbs_thash(n, seed, _hbs_adrs(_HBS_FTS_ROOTS, loc...), roots...)

# The public key of the key at `loc`.
function _fors_pk(n::Int, a::Int, k::Int, sk_seed::Vector{UInt8}, seed::Vector{UInt8}, loc)
    return _fors_compress(n, seed, loc, [merkle_root(_fors_tree(n, a, sk_seed, seed, loc, i)) for i in 0:k - 1])
end

# Sign the digest `md`; the signature has `k * (a + 1) * n` bytes.
function _fors_sign(n::Int, a::Int, k::Int, md::Vector{UInt8}, sk_seed::Vector{UInt8}, seed::Vector{UInt8}, loc)
    idx = _hbs_indices(md, k, a)
    sig = UInt8[]
    for i in 0:k - 1
        levels = _fors_tree(n, a, sk_seed, seed, loc, i)
        append!(sig, _fors_sk(n, sk_seed, loc, i * 2^a + idx[i + 1]))
        append!(sig, reduce(vcat, merkle_auth_path(levels, idx[i + 1])))
    end
    return sig
end

# The public key that a signature of the digest `md` belongs to.
function _fors_pk_from_sig(n::Int, a::Int, k::Int, sig::Vector{UInt8}, md::Vector{UInt8}, seed::Vector{UInt8}, loc)
    @req length(sig) == k * (a + 1) * n "signature must have length $(k * (a + 1) * n), got $(length(sig))"
    idx = _hbs_indices(md, k, a)
    roots = Vector{Vector{UInt8}}()
    for i in 0:k - 1
        parts = _hbs_split(sig[i * (a + 1) * n + 1:(i + 1) * (a + 1) * n], n)
        leaf = _fors_leaf(n, seed, loc, i * 2^a + idx[i + 1], parts[1])
        push!(roots, _merkle_root_from_path(leaf, idx[i + 1], parts[2:end]; node = _fors_node(n, seed, loc, i)))
    end
    return _fors_compress(n, seed, loc, roots)
end
