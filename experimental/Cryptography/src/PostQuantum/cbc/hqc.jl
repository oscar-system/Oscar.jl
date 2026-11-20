# ---------------------------------------------------------------------------
# The public-key encryption scheme of HQC (Hamming Quasi-Cyclic)
# [ABDGZ18](@cite). All computations take place in the ring
# GF(2)[x]/(x^n - 1). The secret key consists of two sparse polynomials x, y,
# the public key of a random polynomial h and s = x + h*y. A ciphertext
# (u, v) hides the codeword of the message under noise that is small enough
# to be removed with the secret key and the concatenated code of HQC.
# ---------------------------------------------------------------------------

struct HQC
    n::Int
    n1::Int
    n2::Int
    k::Int
    w::Int
    wr::Int
    we::Int
    name::String
end

function Base.show(io::IO, p::HQC)
    print(io, "HQC parameters ", p.name, " (n = ", p.n, ", n1 = ", p.n1, ", n2 = ", p.n2, ", k = ", p.k)
    print(io, ", w = ", p.w, ", wr = ", p.wr, ", we = ", p.we, ")")
end

@doc raw"""
    instantiate_hqc(n::Int, n1::Int, k::Int, mult::Int, w::Int, wr::Int, we::Int, name::String) -> HQC

Create a parameter set for the encryption scheme HQC working in the ring
$\mathrm{GF}(2)[x]/(x^n - 1)$. A message block of `k` bytes is encoded with a
Reed-Solomon code over $\mathrm{GF}(2^8)$ of length `n1` and dimension `k`,
whose symbols are encoded with the Reed-Muller code $\mathrm{RM}(1,7)$ of length
128, repeated `mult` times, to a codeword of $n_1 n_2$ bits with
$n_2 = 128 \cdot mult$. The secret polynomials have weight `w`, the random
polynomials of the encryption weights `wr` and `we`. The parameters must
satisfy the following restrictions:
- `n` is prime and `n >= n1 * n2`;
- `1 <= k < n1 <= 255` and `mult >= 1`;
- `w`, `wr` and `we` are between 1 and `n`.

The label `name` is only used for printing. The predefined sets `hqc_toy`,
`hqc_1`, `hqc_3` and `hqc_5` are created this way.

# Examples
```jldoctest
julia> instantiate_hqc(1283, 10, 4, 1, 8, 8, 8, "my-hqc")
HQC parameters my-hqc (n = 1283, n1 = 10, n2 = 128, k = 4, w = 8, wr = 8, we = 8)

julia> instantiate_hqc(1279, 10, 4, 1, 8, 8, 8, "bad")
ERROR: ArgumentError: n must be at least n1 * n2 = 1280, got 1279
[...]
```
"""
function instantiate_hqc(n::Int, n1::Int, k::Int, mult::Int, w::Int, wr::Int, we::Int, name::String)
    @req 1 <= k < n1 <= 255 "k and n1 must satisfy 1 <= k < n1 <= 255, got k = $k and n1 = $n1"
    @req mult >= 1 "mult must be positive, got $mult"
    n2 = 128 * mult
    @req n >= n1 * n2 "n must be at least n1 * n2 = $(n1 * n2), got $n"
    @req is_prime(n) "n must be prime, got $n"
    @req all(v -> 1 <= v <= n, (w, wr, we)) "w, wr and we must be between 1 and n, got $w, $wr and $we"
    return HQC(n, n1, n2, k, w, wr, we, name)
end

# Parameter sets. `hqc_1`, `hqc_3` and `hqc_5` use the sizes of the sets
# HQC-1, HQC-3 and HQC-5 of the NIST submission; the scheme is not compatible
# with it. `hqc_toy` is tiny and not secure.
const hqc_toy = instantiate_hqc(1283, 10, 4, 1, 8, 8, 8, "hqc-toy")
const hqc_1 = instantiate_hqc(17669, 46, 16, 3, 66, 75, 75, "hqc-1")
const hqc_3 = instantiate_hqc(35851, 56, 24, 5, 100, 114, 114, "hqc-3")
const hqc_5 = instantiate_hqc(57637, 90, 32, 5, 131, 149, 149, "hqc-5")

const hqc_all = [hqc_toy, hqc_1, hqc_3, hqc_5]

# Polynomial ring over GF(2) and the modulus x^n - 1 of the parameters `p`.
function _hqc_ring(p::HQC)
    R, x = polynomial_ring(GF(2), :x)
    return R, x^p.n - 1
end

# Polynomial with the coefficients `bits`, starting with the constant term.
_hqc_poly(R::FqPolyRing, bits::AbstractVector{Int}) = R(base_ring(R).(bits))

# Coefficients of the terms of degree less than `len` of the polynomial `f`.
_hqc_bits(f::FqPolyRingElem, len::Int) = [is_one(coeff(f, i)) ? 1 : 0 for i in 0:len-1]

# Random polynomial of degree less than `n` with exactly `w` nonzero terms.
function _hqc_sparse(rng::AbstractRNG, R::FqPolyRing, n::Int, w::Int)
    bits = zeros(Int, n)
    bits[randperm(rng, n)[1:w]] .= 1
    return _hqc_poly(R, bits)
end

# Codeword of the concatenated code for the `k` bytes `msg`, as a polynomial.
function _hqc_encode(p::HQC, R::FqPolyRing, msg::Vector{UInt8})
    mult = p.n2 ÷ 128
    bits = reduce(vcat, [_rm_encode(b, mult) for b in _rs_encode(msg, p.n1)])
    return _hqc_poly(R, bits)
end

# Message of `k` bytes of the noisy codeword `c` of the concatenated code, or
# `nothing` if the Reed-Solomon decoder fails.
function _hqc_decode(p::HQC, c::FqPolyRingElem)
    mult = p.n2 ÷ 128
    bits = _hqc_bits(c, p.n1 * p.n2)
    symbols = [_rm_decode(bits[(i - 1) * p.n2 + 1:i * p.n2], mult) for i in 1:p.n1]
    return _rs_decode(symbols, p.k)
end

@doc raw"""
    keygen(p::HQC; rng::AbstractRNG = RandomDevice()) -> Tuple{Tuple{FqPolyRingElem, FqPolyRingElem}, Tuple{FqPolyRingElem, FqPolyRingElem}}

Generate a key pair for the HQC parameters `p`. The secret key is a pair
`(x, y)` of random polynomials of weight $w$, the public key is the pair
`(h, s)` of a uniformly random polynomial $h$ and $s = x + h y$ in
$\mathrm{GF}(2)[x]/(x^n - 1)$. Return `(pk, sk)`.

# Examples
```jldoctest
julia> (h, s), (x, y) = keygen(hqc_toy);

julia> degree(h) < hqc_toy.n, count(!is_zero, coefficients(x))
(true, 8)
```
"""
function keygen(p::HQC; rng::AbstractRNG = RandomDevice())
    R, f = _hqc_ring(p)
    h = _hqc_poly(R, rand(rng, 0:1, p.n))
    x = _hqc_sparse(rng, R, p.n, p.w)
    y = _hqc_sparse(rng, R, p.n, p.w)
    s = mod(x + h * y, f)
    return (h, s), (x, y)
end

@doc raw"""
    encrypt(p::HQC, pk, msg::AbstractVector{UInt8}; rng::AbstractRNG = RandomDevice()) -> Vector{Tuple{FqPolyRingElem, FqPolyRingElem}}

Encrypt the bytes `msg` under the public key `pk = (h, s)`. The bits of the
message are padded with a one and zeros and split into blocks of $k$ bytes.
Each block $m$ is encoded to a codeword $c$ of the concatenated code and
encrypted with random polynomials $r_1, r_2$ of weight $w_r$ and $e$ of weight
$w_e$ as the pair $u = r_1 + h r_2$ and $v = c + s r_2 + e$, where $v$ is
truncated to $n_1 n_2$ coefficients. The ciphertext is the vector of these
pairs.

# Examples
```jldoctest
julia> pk, sk = keygen(hqc_toy);

julia> ct = encrypt(hqc_toy, pk, UInt8[1, 2, 3, 4]);

julia> length(ct)
2
```
"""
function encrypt(p::HQC, pk, msg::AbstractVector{UInt8}; rng::AbstractRNG = RandomDevice())
    h, s = pk
    R, f = _hqc_ring(p)
    bits = _cbc_pad(_bytes_to_bits(msg), 8 * p.k)
    ct = Tuple{FqPolyRingElem, FqPolyRingElem}[]
    for i in 1:length(bits) ÷ (8 * p.k)
        c = _hqc_encode(p, R, _bits_to_bytes(bits[(i - 1) * 8 * p.k + 1:i * 8 * p.k]))
        r1 = _hqc_sparse(rng, R, p.n, p.wr)
        r2 = _hqc_sparse(rng, R, p.n, p.wr)
        e = _hqc_sparse(rng, R, p.n, p.we)
        u = mod(r1 + h * r2, f)
        v = truncate(c + mod(s * r2, f) + e, p.n1 * p.n2)
        push!(ct, (u, v))
    end
    return ct
end

@doc raw"""
    decrypt(p::HQC, sk, ct) -> Union{Vector{UInt8}, Nothing}

Decrypt the ciphertext `ct` with the secret key `sk = (x, y)`. For each pair
$(u, v)$ the polynomial $v - u y = c + x r_2 - r_1 y + e$, truncated to
$n_1 n_2$ coefficients, is the codeword $c$ of the message block with an error
of low weight, which is removed by decoding the Reed-Muller and then the
Reed-Solomon code. Returns `nothing` if decoding fails or the result is not a
correctly padded message. For the predefined parameter sets, decoding fails
only with negligible probability.

# Examples
```jldoctest
julia> pk, sk = keygen(hqc_toy);

julia> ct = encrypt(hqc_toy, pk, Vector{UInt8}("oscar"));

julia> String(decrypt(hqc_toy, sk, ct))
"oscar"
```
"""
function decrypt(p::HQC, sk, ct)
    x, y = sk
    _, f = _hqc_ring(p)
    bits = Int[]
    for (u, v) in ct
        m = _hqc_decode(p, v - truncate(mod(u * y, f), p.n1 * p.n2))
        m === nothing && return nothing
        append!(bits, _bytes_to_bits(m))
    end
    bits = _cbc_unpad(bits)
    return bits === nothing ? nothing : _bits_to_bytes(bits)
end
