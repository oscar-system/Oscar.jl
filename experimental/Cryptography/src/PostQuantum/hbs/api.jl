# ---------------------------------------------------------------------------
# The interface of the hash-based signature schemes. It matches the one of the
# UOV signature scheme: keys, signatures and messages are byte vectors.
# ---------------------------------------------------------------------------

@doc raw"""
    keygen(p::HBS; rng::AbstractRNG = RandomDevice()) -> Tuple{Vector{UInt8}, Vector{UInt8}}

Generate a key pair for the parameter set `p`, returning the public key `pk`
first and the secret key `sk` second. Both are byte vectors of the lengths
`p.pk_sz` and `p.sk_sz`. For a scheme with a single Merkle tree this computes
the public keys of all $2^h$ Winternitz keys, for the schemes with a hypertree
only those of the single tree in the top layer.

# Examples
```jldoctest
julia> pk, sk = keygen(hbs_mss_toy);

julia> length(pk), length(sk)
(16, 40)

julia> pk, sk = keygen(hbs_fors_toy);

julia> length(pk), length(sk)
(16, 32)
```
"""
function keygen(p::HBS; rng::AbstractRNG = RandomDevice())
    return p.scheme == :mss ? _mss_keygen(p, rng) : _ht_keygen(p, rng)
end

@doc raw"""
    sign(p::HBS, msg::Vector{UInt8}, sk::Vector{UInt8}) -> Vector{UInt8}

Sign the message `msg` with the secret key `sk` of the parameter set `p`,
returning the signature, a byte vector of length `p.sig_sz`.

For a scheme with a single Merkle tree (`p.scheme == :mss`) the scheme is
stateful: the counter at the end of `sk` is increased in place, and at most
$2^h$ messages can be signed with a key. Signing then needs to rebuild the
whole tree, which makes it as expensive as the key generation. The schemes with
a hypertree are stateless, but the few-time keys at the bottom lose security
if the number of signatures comes near to $2^h$.

# Examples
```jldoctest
julia> pk, sk = keygen(hbs_mss_toy);

julia> sig = sign(hbs_mss_toy, UInt8[1, 2, 3], sk);

julia> length(sig) == hbs_mss_toy.sig_sz
true

julia> verify(hbs_mss_toy, sig, UInt8[1, 2, 3], pk)
true
```
"""
function Base.sign(p::HBS, msg::Vector{UInt8}, sk::Vector{UInt8})
    return p.scheme == :mss ? _mss_sign(p, msg, sk) : _ht_sign(p, msg, sk)
end

@doc raw"""
    verify(p::HBS, sig::Vector{UInt8}, msg::Vector{UInt8}, pk::Vector{UInt8}) -> Bool

Return `true` if `sig` is a valid signature of the message `msg` under the
public key `pk` of the parameter set `p`.

# Examples
```jldoctest
julia> pk, sk = keygen(hbs_fors_toy);

julia> msg = Vector{UInt8}("oscar");

julia> sig = sign(hbs_fors_toy, msg, sk);

julia> verify(hbs_fors_toy, sig, msg, pk)
true

julia> verify(hbs_fors_toy, sig, Vector{UInt8}("Oscar"), pk)
false
```
"""
function verify(p::HBS, sig::Vector{UInt8}, msg::Vector{UInt8}, pk::Vector{UInt8})
    return p.scheme == :mss ? _mss_verify(p, sig, msg, pk) : _ht_verify(p, sig, msg, pk)
end

@doc raw"""
    open(p::HBS, sm::Vector{UInt8}, pk::Vector{UInt8})

Recover the message from a signed message `sm` (the message followed by its
signature), verifying it against the public key `pk`. Returns `nothing` if the
signature is invalid.

# Examples
```jldoctest
julia> pk, sk = keygen(hbs_fors_toy);

julia> msg = UInt8[1, 2, 3];

julia> sm = vcat(msg, sign(hbs_fors_toy, msg, sk));

julia> open(hbs_fors_toy, sm, pk)
3-element Vector{UInt8}:
 0x01
 0x02
 0x03
```
"""
function Base.open(p::HBS, sm::Vector{UInt8}, pk::Vector{UInt8})
    length(sm) >= p.sig_sz || return nothing
    msg_sz = length(sm) - p.sig_sz
    msg = sm[1:msg_sz]
    sig = sm[msg_sz + 1:end]
    return verify(p, sig, msg, pk) ? msg : nothing
end
