# ---------------------------------------------------------------------------
# This file implements the SHAKE256 extendable output function, which is part
# of the SHA-3 standard, see [NIS15](@cite). It is implemented as a sponge
# construction over the Keccak-f[1600] permutation with a rate of 136 bytes
# and the "10*1" padding scheme.
# ---------------------------------------------------------------------------

# Round constants of the Keccak-f[1600] permutation, see Table 3 of
# [NIS15](@cite).
const _KECCAK_RC = (
    0x0000000000000001, 0x0000000000008082, 0x800000000000808A, 0x8000000080008000,
    0x000000000000808B, 0x0000000080000001, 0x8000000080008081, 0x8000000000008009,
    0x000000000000008A, 0x0000000000000088, 0x0000000080008009, 0x000000008000000A,
    0x000000008000808B, 0x800000000000008B, 0x8000000000008089, 0x8000000000008003,
    0x8000000000008002, 0x8000000000000080, 0x000000000000800A, 0x800000008000000A,
    0x8000000080008081, 0x8000000000008080, 0x0000000080000001, 0x8000000080008008,
)

# Rotation offsets of the rho step for the 25 lanes of the state, see Table 2
# of [NIS15](@cite). Lane `i` (1-based) is the lane with first coordinate
# `(i-1) mod 5` and second coordinate `(i-1) div 5` of the state.
const _KECCAK_ROT = (0, 1, 62, 28, 27, 36, 44, 6, 55, 20, 3, 10, 43, 25, 39, 41, 45, 15, 21, 8, 18, 2, 61, 56, 14)

# Destination lane (1-based) of the pi permutation for each source lane
# (1-based), as defined in [NIS15](@cite).
const _KECCAK_PI = (1, 11, 21, 6, 16, 17, 2, 12, 22, 7, 8, 18, 3, 13, 23, 24, 9, 19, 4, 14, 15, 25, 10, 20, 5)

# Apply the Keccak-f[1600] permutation to the 25-lane state `s` in place,
# using `b` as a 25-lane scratch buffer.
function _keccak_f1600!(s::Vector{UInt64}, b::Vector{UInt64})
    @inbounds for round in 1:24
        # theta
        c1 = s[1] ⊻ s[6] ⊻ s[11] ⊻ s[16] ⊻ s[21]
        c2 = s[2] ⊻ s[7] ⊻ s[12] ⊻ s[17] ⊻ s[22]
        c3 = s[3] ⊻ s[8] ⊻ s[13] ⊻ s[18] ⊻ s[23]
        c4 = s[4] ⊻ s[9] ⊻ s[14] ⊻ s[19] ⊻ s[24]
        c5 = s[5] ⊻ s[10] ⊻ s[15] ⊻ s[20] ⊻ s[25]
        d1 = c5 ⊻ bitrotate(c2, 1)
        d2 = c1 ⊻ bitrotate(c3, 1)
        d3 = c2 ⊻ bitrotate(c4, 1)
        d4 = c3 ⊻ bitrotate(c5, 1)
        d5 = c4 ⊻ bitrotate(c1, 1)
        s[1] ⊻= d1; s[6] ⊻= d1; s[11] ⊻= d1; s[16] ⊻= d1; s[21] ⊻= d1
        s[2] ⊻= d2; s[7] ⊻= d2; s[12] ⊻= d2; s[17] ⊻= d2; s[22] ⊻= d2
        s[3] ⊻= d3; s[8] ⊻= d3; s[13] ⊻= d3; s[18] ⊻= d3; s[23] ⊻= d3
        s[4] ⊻= d4; s[9] ⊻= d4; s[14] ⊻= d4; s[19] ⊻= d4; s[24] ⊻= d4
        s[5] ⊻= d5; s[10] ⊻= d5; s[15] ⊻= d5; s[20] ⊻= d5; s[25] ⊻= d5

        # rho & pi
        for i in 1:25
            b[_KECCAK_PI[i]] = bitrotate(s[i], _KECCAK_ROT[i])
        end

        # chi
        for base in (1, 6, 11, 16, 21)
            s[base] = b[base] ⊻ (~b[base+1] & b[base+2])
            s[base+1] = b[base+1] ⊻ (~b[base+2] & b[base+3])
            s[base+2] = b[base+2] ⊻ (~b[base+3] & b[base+4])
            s[base+3] = b[base+3] ⊻ (~b[base+4] & b[base])
            s[base+4] = b[base+4] ⊻ (~b[base] & b[base+1])
        end

        # iota
        s[1] ⊻= _KECCAK_RC[round]
    end
    return s
end

# Read 8 consecutive bytes of `msg` starting at position `p` as a
# little-endian 64-bit lane.
@inline function _keccak_read_lane(msg::AbstractVector{UInt8}, p::Int)::UInt64
    @inbounds return UInt64(msg[p]) | (UInt64(msg[p+1]) << 8) | (UInt64(msg[p+2]) << 16) |
                       (UInt64(msg[p+3]) << 24) | (UInt64(msg[p+4]) << 32) |
                       (UInt64(msg[p+5]) << 40) | (UInt64(msg[p+6]) << 48) |
                       (UInt64(msg[p+7]) << 56)
end

@doc raw"""
    shake_256(msg::AbstractVector{UInt8}, outlen::Int) -> Vector{UInt8}

Return the first `outlen` bytes of the SHAKE256 extendable output function
applied to the byte vector `msg`, see [NIS15](@cite) for details.

# Examples
```jldoctest
julia> bytes2hex(shake_256(UInt8[], 32))
"46b9dd2b0ba88d13233b3feb743eeb243fcd52ea62b81b82b50c27646ed5762f"
```
"""
function shake_256(msg::AbstractVector{UInt8}, outlen::Int)::Vector{UInt8}
    @req outlen >= 0 "output length must be non-negative, got $outlen"

    rate = 136
    s = zeros(UInt64, 25)
    b = Vector{UInt64}(undef, 25)

    # absorb complete rate blocks
    n = length(msg)
    p = 1
    while p + rate - 1 <= n
        @inbounds for i in 1:17
            s[i] ⊻= _keccak_read_lane(msg, p + 8(i-1))
        end
        _keccak_f1600!(s, b)
        p += rate
    end

    # absorb the remaining bytes (0 <= r <= rate - 1)
    r = n - p + 1
    nl = r ÷ 8
    @inbounds for i in 1:nl
        s[i] ⊻= _keccak_read_lane(msg, p + 8(i-1))
    end
    rb = r - 8nl
    if rb > 0
        v = UInt64(0)
        @inbounds for j in 1:rb
            v |= UInt64(msg[p + 8nl + j - 1]) << (8(j-1))
        end
        s[nl+1] ⊻= v
    end

    # apply the pad10*1 padding with the SHAKE domain separation: the first
    # padding byte 0x1f goes right after the message, and the final bit 0x80
    # into the last position of the block
    off = r
    s[off ÷ 8 + 1] ⊻= UInt64(0x1f) << (8 * (off % 8))
    s[17] ⊻= 0x8000000000000000
    _keccak_f1600!(s, b)

    # squeeze
    out = Vector{UInt8}(undef, outlen)
    o = 1
    while o <= outlen
        take = min(rate, outlen - o + 1)
        nl = take ÷ 8
        @inbounds for i in 1:nl
            lane = s[i]
            for j in 1:8
                out[o + 8(i-1) + j - 1] = UInt8((lane >> (8(j-1))) & 0xff)
            end
        end
        rb = take - 8nl
        if rb > 0
            lane = s[nl + 1]
            @inbounds for j in 1:rb
                out[o + 8nl + j - 1] = UInt8((lane >> (8(j-1))) & 0xff)
            end
        end
        o += take
        if o <= outlen
            _keccak_f1600!(s, b)
        end
    end
    return out
end
