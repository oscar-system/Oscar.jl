# ---------------------------------------------------------------------------
# Lamport's one-time signature scheme [Lam79](@cite). It is the simplest hash-based
# signature and mainly meant for teaching: a key pair must only be used to sign
# one message. The message is hashed to `8n` bits with SHAKE256, and each bit
# is signed by revealing one of two secret strings of `n` bytes.
# ---------------------------------------------------------------------------

# Hash of a message to the bits that are signed, as a vector of 0/1 values.
_lamport_bits(n::Int, msg::Vector{UInt8}) = [Int((b >> (7 - i)) & 1) for b in shake_256(msg, n) for i in 0:7]

@doc raw"""
    lamport_keygen(n::Int = 16; rng::AbstractRNG = RandomDevice()) -> Tuple{Vector{UInt8}, Vector{UInt8}}

Generate a key pair of Lamport's one-time signature scheme with hash length `n`
bytes. The secret key consists of $2 \cdot 8n$ random strings of `n` bytes, the
public key consists of their hashes, both stored as one byte vector. Return
`(pk, sk)`. A key pair must be used for one signature only.

# Examples
```jldoctest
julia> pk, sk = lamport_keygen(4);

julia> length(pk), length(sk)
(256, 256)
```
"""
function lamport_keygen(n::Int = 16; rng::AbstractRNG = RandomDevice())
    @req n > 0 "n must be positive, got $n"
    sk = rand(rng, UInt8, 16n * n)
    pk = reduce(vcat, [shake_256(sk[(i - 1) * n + 1:i * n], n) for i in 1:16n])
    return pk, sk
end

@doc raw"""
    lamport_sign(n::Int, msg::Vector{UInt8}, sk::Vector{UInt8}) -> Vector{UInt8}

Sign `msg` with the secret key `sk`. For the $i$-th bit $b_i$ of the hash of
`msg` the string number $2i + b_i$ of the secret key is revealed, so the
signature has $8n \cdot n$ bytes.

# Examples
```jldoctest
julia> pk, sk = lamport_keygen(4);

julia> length(lamport_sign(4, UInt8[1, 2, 3], sk))
128
```
"""
function lamport_sign(n::Int, msg::Vector{UInt8}, sk::Vector{UInt8})
    @req length(sk) == 16n * n "secret key must have length $(16n * n), got $(length(sk))"
    bits = _lamport_bits(n, msg)
    return reduce(vcat, [sk[(2(i - 1) + b) * n + 1:(2(i - 1) + b + 1) * n] for (i, b) in enumerate(bits)])
end

@doc raw"""
    lamport_verify(n::Int, sig::Vector{UInt8}, msg::Vector{UInt8}, pk::Vector{UInt8}) -> Bool

Return `true` if `sig` is a valid Lamport signature of `msg` for the public key
`pk`, that is, if the hash of each revealed string is the matching entry of `pk`.

# Examples
```jldoctest
julia> pk, sk = lamport_keygen(4);

julia> sig = lamport_sign(4, UInt8[1, 2, 3], sk);

julia> lamport_verify(4, sig, UInt8[1, 2, 3], pk)
true

julia> lamport_verify(4, sig, UInt8[1, 2, 4], pk)
false
```
"""
function lamport_verify(n::Int, sig::Vector{UInt8}, msg::Vector{UInt8}, pk::Vector{UInt8})
    @req length(sig) == 8n * n "signature must have length $(8n * n), got $(length(sig))"
    @req length(pk) == 16n * n "public key must have length $(16n * n), got $(length(pk))"
    bits = _lamport_bits(n, msg)
    for (i, b) in enumerate(bits)
        j = 2(i - 1) + b
        shake_256(sig[(i - 1) * n + 1:i * n], n) == pk[j * n + 1:(j + 1) * n] || return false
    end
    return true
end
