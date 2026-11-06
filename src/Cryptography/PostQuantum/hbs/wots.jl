# ---------------------------------------------------------------------------
# The Winternitz one-time signature scheme in the variant WOTS+ [Hul13](@cite).
# A message of `n` bytes is written in base `w`, and a checksum is appended.
# Each digit is signed by the matching number of steps in a hash chain that
# starts at a secret value. The public key is the compression of the chain ends.
#
# A key is identified by `loc = (layer, tree, keypair)`, its position in the
# scheme. All secret values are derived from `sk_seed`, all hashing is bound
# to the public seed `seed`.
# ---------------------------------------------------------------------------

# Number of chains `(len1, len2)` for the message and the checksum.
function _wots_lens(n::Int, w::Int)
    lw = trailing_zeros(w)
    len1 = cld(8n, lw)
    len2 = ndigits(len1 * (w - 1); base = w)
    return len1, len2
end

_wots_len(n::Int, w::Int) = sum(_wots_lens(n, w))

# The base `w` digits of the message `m` followed by the digits of the checksum.
function _wots_digits(n::Int, w::Int, m::Vector{UInt8})
    @req length(m) == n "message must have length $n, got $(length(m))"
    lw = trailing_zeros(w)
    len1, len2 = _wots_lens(n, w)
    bits = vcat([Int((b >> (7 - i)) & 1) for b in m for i in 0:7], zeros(Int, len1 * lw - 8n))
    d = [evalpoly(2, reverse(bits[(j - 1) * lw + 1:j * lw])) for j in 1:len1]
    csum = sum(w - 1 - x for x in d)
    return vcat(d, reverse(digits(csum; base = w, pad = len2)))
end

# Apply `steps` steps of the chain number `chain` to `x`, starting at step `start`.
function _wots_chain(n::Int, seed::Vector{UInt8}, loc, chain::Int, x::Vector{UInt8}, start::Int, steps::Int)
    for j in start:start + steps - 1
        x = _hbs_thash(n, seed, _hbs_adrs(_HBS_WOTS_HASH, loc..., chain, j), x)
    end
    return x
end

# The secret value at the beginning of the chain number `chain`.
_wots_sk(n::Int, sk_seed::Vector{UInt8}, loc, chain::Int) =
    _hbs_prf(n, sk_seed, _hbs_adrs(_HBS_PRF, _HBS_WOTS_HASH, loc..., chain))

# Compress the chain ends to the public key.
_wots_compress(n::Int, seed::Vector{UInt8}, loc, ends::Vector{Vector{UInt8}}) =
    _hbs_thash(n, seed, _hbs_adrs(_HBS_WOTS_PK, loc...), ends...)

# The public key of the key at `loc`.
function _wots_pk(n::Int, w::Int, sk_seed::Vector{UInt8}, seed::Vector{UInt8}, loc)
    ends = [_wots_chain(n, seed, loc, i, _wots_sk(n, sk_seed, loc, i), 0, w - 1) for i in 0:_wots_len(n, w) - 1]
    return _wots_compress(n, seed, loc, ends)
end

# Sign the message `m` of `n` bytes; the signature is a vector of `len` strings of `n` bytes.
function _wots_sign(n::Int, w::Int, m::Vector{UInt8}, sk_seed::Vector{UInt8}, seed::Vector{UInt8}, loc)
    d = _wots_digits(n, w, m)
    return [_wots_chain(n, seed, loc, i - 1, _wots_sk(n, sk_seed, loc, i - 1), 0, d[i]) for i in eachindex(d)]
end

# The public key that a signature of the message `m` belongs to.
function _wots_pk_from_sig(n::Int, w::Int, sig::Vector{Vector{UInt8}}, m::Vector{UInt8}, seed::Vector{UInt8}, loc)
    d = _wots_digits(n, w, m)
    @req length(sig) == length(d) "signature must have $(length(d)) parts, got $(length(sig))"
    ends = [_wots_chain(n, seed, loc, i - 1, sig[i], d[i], w - 1 - d[i]) for i in eachindex(d)]
    return _wots_compress(n, seed, loc, ends)
end
