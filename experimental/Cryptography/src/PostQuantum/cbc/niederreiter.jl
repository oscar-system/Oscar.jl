# ---------------------------------------------------------------------------
# Niederreiter's public-key encryption scheme [Nie86](@cite) with binary Goppa
# codes, the dual version of McEliece's scheme. The public key is a scrambled
# parity check matrix S*H*P of a random Goppa code; a message block is encoded
# as an error vector of weight t, whose syndrome is the ciphertext. Only the
# owner of the code can find the error vector from its syndrome.
# ---------------------------------------------------------------------------

struct Niederreiter
    m::Int
    n::Int
    t::Int
    l::Int
    name::String
end

function Base.show(io::IO, p::Niederreiter)
    print(io, "Niederreiter parameters ", p.name, " (m = ", p.m, ", n = ", p.n, ", t = ", p.t, ", l = ", p.l, ")")
end

@doc raw"""
    instantiate_niederreiter(m::Int, n::Int, t::Int, name::String) -> Niederreiter

Create a parameter set for Niederreiter's encryption scheme with binary Goppa
codes of length `n` over $\mathrm{GF}(2^m)$ that correct `t` errors. A block
of a message has $l = \lfloor \log_2 \binom{n}{t} \rfloor$ bits, so that it can
be encoded as a binary vector of length $n$ and weight $t$. The parameters must
satisfy the following restrictions:
- `2 <= m <= 16`;
- `t >= 2`;
- `m * t < n <= 2^m`.

The label `name` is only used for printing. The predefined sets
`niederreiter_toy`, `niederreiter_light` and `niederreiter_standard` are
created this way.

# Examples
```jldoctest
julia> instantiate_niederreiter(6, 64, 5, "my-niederreiter")
Niederreiter parameters my-niederreiter (m = 6, n = 64, t = 5, l = 22)

julia> instantiate_niederreiter(6, 64, 1, "bad")
ERROR: ArgumentError: t must be at least 2, got 1
[...]
```
"""
function instantiate_niederreiter(m::Int, n::Int, t::Int, name::String)
    _check_goppa_parameters(m, n, t)
    l = ndigits(binomial(big(n), t); base = 2) - 1
    return Niederreiter(m, n, t, l, name)
end

# Parameter sets with the same codes as the ones of McEliece's scheme.
const niederreiter_toy = instantiate_niederreiter(6, 64, 5, "niederreiter-toy")
const niederreiter_light = instantiate_niederreiter(10, 1024, 50, "niederreiter-light")
const niederreiter_standard = instantiate_niederreiter(12, 3488, 64, "niederreiter-standard")

const niederreiter_all = [niederreiter_toy, niederreiter_light, niederreiter_standard]

# Positions (between 1 and n) of the ones of the binary vector of length `n`
# and weight `t` with number `N`, where 0 <= N < binomial(n, t). It uses the
# combinatorial number system: N is written uniquely as
# binomial(c_t, t) + ... + binomial(c_1, 1) with n > c_t > ... > c_1 >= 0.
function _constant_weight_encode(N::BigInt, n::Int, t::Int)
    pos = Int[]
    c = n
    for i in t:-1:1
        c -= 1
        b = binomial(big(c), i)
        while b > N
            b = b * (c - i) ÷ c
            c -= 1
        end
        push!(pos, c + 1)
        N -= b
    end
    return reverse(pos)
end

# Number of the binary vector with the ones at the sorted positions `pos`,
# the inverse of `_constant_weight_encode`.
function _constant_weight_decode(pos::Vector{Int})
    return sum(binomial(big(c - 1), i) for (i, c) in enumerate(pos); init = big(0))
end

@doc raw"""
    keygen(p::Niederreiter; rng::AbstractRNG = RandomDevice()) -> Tuple{FqMatrix, Tuple{BinaryGoppaCode, FqMatrix, Vector{Int}}}

Generate a key pair for the Niederreiter parameters `p`. A random binary Goppa
code with binary parity check matrix $H$ is chosen together with a random
invertible $mt \times mt$ matrix $S$ and a random permutation matrix $P$ over
$\mathrm{GF}(2)$. The public key is the $mt \times n$ matrix $SHP$, the secret
key is the triple `(C, S^{-1}, perm)` of the code, the inverse of $S$ and the
permutation. Return `(pk, sk)`.

# Examples
```jldoctest
julia> pk, sk = keygen(niederreiter_toy);

julia> size(pk)
(30, 64)
```
"""
function keygen(p::Niederreiter; rng::AbstractRNG = RandomDevice())
    C = _goppa_code(rng, p.m, p.n, p.t)
    S = _random_invertible_gf2(rng, p.m * p.t)
    perm = randperm(rng, p.n)
    pk = (S * C.H)[:, perm]
    return pk, (C, inv(S), perm)
end

@doc raw"""
    encrypt(p::Niederreiter, pk::FqMatrix, msg::AbstractVector{UInt8}) -> FqMatrix

Encrypt the bytes `msg` under the public key `pk`. The bits of the message are
padded with a one and zeros and split into blocks of $l$ bits. Each block is
read as a number below $\binom{n}{t}$ and encoded as a binary vector $e$ of
length $n$ and weight $t$, which is encrypted as its syndrome $pk \cdot e^T$.
The ciphertext is the matrix over $\mathrm{GF}(2)$ with these syndromes as rows.
The encryption is deterministic; the keyword `rng` is accepted for
compatibility with the other schemes and not used.

# Examples
```jldoctest
julia> pk, sk = keygen(niederreiter_toy);

julia> ct = encrypt(niederreiter_toy, pk, UInt8[1, 2, 3]);

julia> size(ct)
(2, 30)
```
"""
function encrypt(p::Niederreiter, pk::FqMatrix, msg::AbstractVector{UInt8}; rng::AbstractRNG = RandomDevice())
    bits = _cbc_pad(_bytes_to_bits(msg), p.l)
    b = length(bits) ÷ p.l
    E = zeros(Int, b, p.n)
    for i in 1:b
        block = bits[(i - 1) * p.l + 1:i * p.l]
        N = sum(big(2)^(j - 1) * block[j] for j in 1:p.l)
        E[i, _constant_weight_encode(N, p.n, p.t)] .= 1
    end
    return matrix(GF(2), E) * transpose(pk)
end

@doc raw"""
    decrypt(p::Niederreiter, sk, ct::FqMatrix) -> Union{Vector{UInt8}, Nothing}

Decrypt the ciphertext `ct` with the secret key `sk`. Multiplying a row of
`ct` with $S^{-1}$ gives the syndrome of the permuted error vector with
respect to the parity check matrix of the Goppa code, from which Patterson's
algorithm [Pat75](@cite) recovers the error vector and thus the message block.
Returns `nothing` if a row cannot be decoded or the result is not a correctly
padded message.

# Examples
```jldoctest
julia> pk, sk = keygen(niederreiter_toy);

julia> ct = encrypt(niederreiter_toy, pk, Vector{UInt8}("oscar"));

julia> String(decrypt(niederreiter_toy, sk, ct))
"oscar"
```
"""
function decrypt(p::Niederreiter, sk, ct::FqMatrix)
    C, S_inv, perm = sk
    @req ncols(ct) == p.m * p.t "ciphertext does not match the parameters"
    syndromes = _gf2_int(S_inv * transpose(ct))
    iperm = invperm(perm)
    bits = Int[]
    for i in 1:nrows(ct)
        pos = _goppa_decode(C, syndromes[:, i])
        (pos === nothing || length(pos) != p.t) && return nothing
        N = _constant_weight_decode(sort(iperm[pos]))
        N < big(2)^p.l || return nothing
        append!(bits, digits(N; base = 2, pad = p.l))
    end
    bits = _cbc_unpad(bits)
    return bits === nothing ? nothing : _bits_to_bytes(bits)
end
