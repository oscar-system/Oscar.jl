# ---------------------------------------------------------------------------
# This file implements Regev's public-key encryption scheme based on the
# learning with errors (LWE) problem, see [Reg05](@cite). Messages are
# encrypted bit by bit. All arithmetic is done with integers modulo q.
#
# The implementation is meant for experiments and teaching. It is neither
# constant-time nor protected against chosen-ciphertext attacks.
# ---------------------------------------------------------------------------
using Random: AbstractRNG, RandomDevice

# ---------------------------------------------------------------------------
# LWE parameter structure
# ---------------------------------------------------------------------------

struct LWE
    n::Int
    m::Int
    q::Int
    eta::Int
    name::String
end

function Base.show(io::IO, p::LWE)
    print(io, "LWE parameters ", p.name, " (n = ", p.n, ", m = ", p.m, ", q = ", p.q, ", eta = ", p.eta, ")")
end

@doc raw"""
    instantiate_lwe(n::Int, m::Int, q::Int, eta::Int, name::String) -> LWE

Create a parameter set for Regev's LWE encryption scheme with secret dimension
$n$, $m$ samples, modulus $q$ and noise drawn from the centered binomial
distribution with parameter $\eta$ (with values in $[-\eta, \eta]$). The
parameters must satisfy the following restrictions:
- `n`, `m`, `eta` are positive;
- `q >= 8`;
- the worst-case decryption noise `m * eta` is less than `q / 4`, so that
  decryption is always correct.

The label `name` is only used for printing. The predefined sets `lwe_toy`,
`lwe_light` and `lwe_standard` are created this way.

# Examples
```jldoctest
julia> p = instantiate_lwe(16, 32, 1031, 2, "my-lwe")
LWE parameters my-lwe (n = 16, m = 32, q = 1031, eta = 2)

julia> instantiate_lwe(16, 32, 100, 2, "bad")
ERROR: ArgumentError: worst-case noise m*eta = 64 must be less than q/4 = 25.0
[...]
```
"""
function instantiate_lwe(n::Int, m::Int, q::Int, eta::Int, name::String)
    @req n > 0 && m > 0 && eta > 0 "n, m and eta must be positive, got n=$n, m=$m, eta=$eta"
    @req q >= 8 "q must be at least 8, got $q"
    @req 4 * m * eta < q "worst-case noise m*eta = $(m*eta) must be less than q/4 = $(q/4)"
    return LWE(n, m, q, eta, name)
end

# Parameter sets. These are illustrative and are not claimed to provide any
# concrete security level, except for `lwe_standard` which follows the usual
# Regev-style choice of dimensions.
const lwe_toy = instantiate_lwe(16, 32, 1031, 2, "lwe-toy")
const lwe_light = instantiate_lwe(128, 256, 4099, 2, "lwe-light")
const lwe_standard = instantiate_lwe(512, 1024, 12289, 2, "lwe-standard")

const lwe_all = [lwe_toy, lwe_light, lwe_standard]

# ---------------------------------------------------------------------------
# Sampling and helpers
# ---------------------------------------------------------------------------

# Sample `k` values of a centered binomial distribution with parameter `eta`.
function _cbd(rng::AbstractRNG, eta::Int, k::Int)
    return [count(rand(rng, Bool, eta)) - count(rand(rng, Bool, eta)) for _ in 1:k]
end

# Reduce `x` into the range [0, q).
_modq(x, q::Int) = mod.(x, q)

# Distance of `x` in Z/qZ to the nearest multiple of q.
_dist0(x::Int, q::Int) = min(x, q - x)

# Convert a vector of bytes to a vector of bits (least significant bit first).
_bytes_to_bits(msg::AbstractVector{UInt8}) = [Int((b >> i) & 1) for b in msg for i in 0:7]

# Convert a vector of bits (length divisible by 8) to bytes.
function _bits_to_bytes(bits::AbstractVector{Int})
    return UInt8[sum(bits[8(j-1)+i] << (i-1) for i in 1:8) for j in 1:length(bits)÷8]
end

# ---------------------------------------------------------------------------
# Key generation, encryption and decryption
# ---------------------------------------------------------------------------

@doc raw"""
    keygen(p::LWE; rng::AbstractRNG = RandomDevice()) -> Tuple{Tuple{Matrix{Int}, Vector{Int}}, Vector{Int}}

Generate a key pair for the LWE parameters `p`. The public key is a pair
`(A, b)` where `A` is a uniformly random $m \times n$ matrix modulo $q$ and
$b = As + e \bmod q$ for a secret vector `s` and a small noise vector `e`.
The secret key is `s`. Return `(pk, sk)`.

# Examples
```jldoctest
julia> (A, b), s = keygen(lwe_toy);

julia> size(A), length(b), length(s)
((32, 16), 32, 16)
```
"""
function keygen(p::LWE; rng::AbstractRNG = RandomDevice())
    A = rand(rng, 0:p.q-1, p.m, p.n)
    s = rand(rng, 0:p.q-1, p.n)
    e = _cbd(rng, p.eta, p.m)
    b = _modq(A * s + e, p.q)
    return (A, b), s
end

@doc raw"""
    encrypt(p::LWE, pk, msg::AbstractVector{UInt8}; rng::AbstractRNG = RandomDevice()) -> Tuple{Matrix{Int}, Vector{Int}}

Encrypt the bytes `msg` under the public key `pk` bit by bit. For each bit
$\mu$ a random $r \in \{0,1\}^m$ is chosen and the pair
$(A^T r, b^T r + \mu \lfloor q/2 \rfloor)$ modulo $q$ is formed. The
ciphertext `(U, v)` consists of the matrix `U` with the vectors $A^T r$ as
columns and the vector `v` of the second components.

# Examples
```jldoctest
julia> pk, sk = keygen(lwe_toy);

julia> U, v = encrypt(lwe_toy, pk, UInt8[1, 2, 3]);

julia> size(U), length(v)
((16, 24), 24)
```
"""
function encrypt(p::LWE, pk, msg::AbstractVector{UInt8}; rng::AbstractRNG = RandomDevice())
    A, b = pk
    bits = _bytes_to_bits(msg)
    R = rand(rng, 0:1, p.m, length(bits))
    U = _modq(transpose(A) * R, p.q)
    v = _modq(vec(transpose(b) * R) .+ (p.q ÷ 2) .* bits, p.q)
    return U, v
end

@doc raw"""
    decrypt(p::LWE, sk, ct) -> Vector{UInt8}

Decrypt the ciphertext `ct = (U, v)` with the secret key `sk`. A bit is
decoded as $1$ if $v_i - s^T U_i \bmod q$ is closer to $\lfloor q/2 \rfloor$
than to $0$. Since the noise is bounded by $m \eta < q/4$ for all parameter
sets created by `instantiate_lwe`, decryption is always correct.

# Examples
```jldoctest
julia> pk, sk = keygen(lwe_toy);

julia> ct = encrypt(lwe_toy, pk, Vector{UInt8}("oscar"));

julia> String(decrypt(lwe_toy, sk, ct))
"oscar"
```
"""
function decrypt(p::LWE, sk, ct)
    U, v = ct
    @req size(U, 1) == p.n && size(U, 2) == length(v) "ciphertext does not match the parameters"
    @req length(v) % 8 == 0 "ciphertext length must be a multiple of 8 bits"
    d = _modq(v - vec(transpose(sk) * U), p.q)
    bits = [Int(_dist0(x, p.q) > p.q ÷ 4) for x in d]
    return _bits_to_bytes(bits)
end
