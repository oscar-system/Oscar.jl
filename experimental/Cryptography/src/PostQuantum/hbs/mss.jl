# ---------------------------------------------------------------------------
# The Merkle signature scheme with a single tree [Mer90](@cite). The public key
# is the root of a tree over the public keys of $2^h$ Winternitz keys. Every
# message is signed with the next unused Winternitz key, so the scheme is
# stateful: the secret key contains a counter, which is increased in place by
# each signature. The same key must never be used twice, so a secret key must
# not be copied and used from two places.
#
# Layout of the keys:
# - public key: `pk_seed || root`;
# - secret key: `sk_seed || sk_prf || pk_seed || root || counter`, where the
#   counter has 8 bytes (little endian);
# - signature: `counter || R || wots signature || authentication path`, where
#   `R` is the randomizer of the message.
# ---------------------------------------------------------------------------

function _mss_keygen(p::HBS, rng::AbstractRNG)
    n = p.n
    sk_seed, sk_prf, pk_seed = rand(rng, UInt8, n), rand(rng, UInt8, n), rand(rng, UInt8, n)
    root = merkle_root(_xmss_tree(n, p.w, p.h, sk_seed, pk_seed, 0, 0))
    return vcat(pk_seed, root), vcat(sk_seed, sk_prf, pk_seed, root, _hbs_adrs(0))
end

# The signature uses the key at the current counter, which is increased in `sk`.
function _mss_sign(p::HBS, msg::Vector{UInt8}, sk::Vector{UInt8})
    n = p.n
    @req length(sk) == p.sk_sz "secret key must have length $(p.sk_sz), got $(length(sk))"
    sk_seed, sk_prf, pk_seed, root = sk[1:n], sk[n + 1:2n], sk[2n + 1:3n], sk[3n + 1:4n]
    counter = sum(Int(sk[4n + i]) << (8(i - 1)) for i in 1:8)
    @req counter < 2^p.h "all $(2^p.h) signatures of this key have been used"
    sk[4n + 1:4n + 8] = _hbs_adrs(counter + 1)
    r = _hbs_prf_msg(n, sk_prf, msg, counter)
    m = _hbs_hmsg(n, r, pk_seed, root, msg)
    levels = _xmss_tree(n, p.w, p.h, sk_seed, pk_seed, 0, 0)
    return vcat(_hbs_adrs(counter), r, _xmss_sign(n, p.w, levels, m, sk_seed, pk_seed, 0, 0, counter))
end

function _mss_verify(p::HBS, sig::Vector{UInt8}, msg::Vector{UInt8}, pk::Vector{UInt8})
    n = p.n
    @req length(sig) == p.sig_sz "signature must have length $(p.sig_sz), got $(length(sig))"
    @req length(pk) == p.pk_sz "public key must have length $(p.pk_sz), got $(length(pk))"
    pk_seed, root = pk[1:n], pk[n + 1:2n]
    counter = sum(Int(sig[i]) << (8(i - 1)) for i in 1:8)
    counter < 2^p.h || return false
    r = sig[9:8 + n]
    m = _hbs_hmsg(n, r, pk_seed, root, msg)
    return _xmss_root_from_sig(n, p.w, p.h, sig[9 + n:end], m, pk_seed, 0, 0, counter) == root
end
