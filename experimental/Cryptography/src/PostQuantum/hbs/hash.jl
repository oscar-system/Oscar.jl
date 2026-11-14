# ---------------------------------------------------------------------------
# Hash functions with domain separation. Every call of the hash function gets
# a public seed and an address, which encodes the type and the position of the
# hash call inside of the scheme. This way no two hash calls of a scheme use
# the same input, which is the usual way to avoid multi-target attacks.
# ---------------------------------------------------------------------------

# Types of addresses.
const _HBS_WOTS_HASH = 0  # chain step of a Winternitz key
const _HBS_WOTS_PK = 1    # compression of the chain ends of a Winternitz key
const _HBS_TREE = 2       # inner node of a Merkle tree over Winternitz keys
const _HBS_FTS_HASH = 3   # leaf of a few-time signature key
const _HBS_FTS_TREE = 4   # inner node of the tree of a few-time signature key
const _HBS_FTS_ROOTS = 5  # compression of the roots of the trees of a FORS key
const _HBS_PRF = 6        # derivation of secret values

# Encode non-negative integers as bytes, using 8 bytes (little endian) each.
function _hbs_adrs(fields::Integer...)
    return UInt8[(UInt64(f) >> (8i)) % UInt8 for f in fields for i in 0:7]
end

# Hash `parts` to `n` bytes, bound to the seed `seed` and the address `adrs`.
function _hbs_thash(n::Int, seed::Vector{UInt8}, adrs::Vector{UInt8}, parts::Vector{UInt8}...)
    return shake_256(vcat(seed, adrs, parts...), n)
end

# Derive the `n` byte secret value at address `adrs` from the secret seed.
function _hbs_prf(n::Int, sk_seed::Vector{UInt8}, adrs::Vector{UInt8})
    return shake_256(vcat(sk_seed, adrs), n)
end

# Randomizer of a message, derived from the secret PRF key.
function _hbs_prf_msg(n::Int, sk_prf::Vector{UInt8}, msg::Vector{UInt8}, counter::Integer = 0)
    return shake_256(vcat(sk_prf, _hbs_adrs(counter), msg), n)
end

# Digest of a message of `outlen` bytes, bound to the randomizer and the public key.
function _hbs_hmsg(outlen::Int, r::Vector{UInt8}, pk_seed::Vector{UInt8}, pk_root::Vector{UInt8}, msg::Vector{UInt8})
    return shake_256(vcat(r, pk_seed, pk_root, msg), outlen)
end

# Read `count` integers of `bits` bits each from `bytes`, most significant bit first.
function _hbs_indices(bytes::Vector{UInt8}, count::Int, bits::Int)
    @req 8 * length(bytes) >= count * bits "not enough bits to read $count integers of $bits bits"
    idx = Vector{Int}(undef, count)
    for j in 1:count
        v = 0
        for b in (j - 1) * bits:j * bits - 1
            v = 2v + Int((bytes[b ÷ 8 + 1] >> (7 - b % 8)) & 1)
        end
        idx[j] = v
    end
    return idx
end

# Split `bytes` into strings of `n` bytes.
function _hbs_split(bytes::Vector{UInt8}, n::Int)
    @req length(bytes) % n == 0 "length $(length(bytes)) is not a multiple of $n"
    return [bytes[(i - 1) * n + 1:i * n] for i in 1:length(bytes) ÷ n]
end
