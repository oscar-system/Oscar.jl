# ---------------------------------------------------------------------------
# McEliece's public-key encryption scheme [McE78](@cite) with binary Goppa
# codes. The public key is a scrambled generator matrix S*G*P of a random
# Goppa code; a message block is encoded with it and t random errors are
# added. Only the owner of the code can remove the errors.
# ---------------------------------------------------------------------------

struct McEliece
    m::Int
    n::Int
    t::Int
    k::Int
    name::String
end

function Base.show(io::IO, p::McEliece)
    print(io, "McEliece parameters ", p.name, " (m = ", p.m, ", n = ", p.n, ", t = ", p.t, ", k = ", p.k, ")")
end

@doc raw"""
    instantiate_mceliece(m::Int, n::Int, t::Int, name::String) -> McEliece

Create a parameter set for McEliece's encryption scheme with binary Goppa
codes of length `n` over $\mathrm{GF}(2^m)$ that correct `t` errors. The codes
have dimension $k = n - mt$. The parameters must satisfy the following
restrictions:
- `2 <= m <= 16`;
- `t >= 2`;
- `m * t < n <= 2^m`.

The label `name` is only used for printing. The predefined sets
`mceliece_toy`, `mceliece_light` and `mceliece_standard` are created this way.

# Examples
```jldoctest
julia> instantiate_mceliece(6, 64, 5, "my-mceliece")
McEliece parameters my-mceliece (m = 6, n = 64, t = 5, k = 34)

julia> instantiate_mceliece(6, 100, 5, "bad")
ERROR: ArgumentError: n must satisfy m*t = 30 < n <= 2^m = 64, got 100
[...]
```
"""
function instantiate_mceliece(m::Int, n::Int, t::Int, name::String)
    _check_goppa_parameters(m, n, t)
    return McEliece(m, n, t, n - m * t, name)
end

# Parameter sets. `mceliece_light` uses the parameters of McEliece's original
# proposal, `mceliece_standard` the code sizes of the set mceliece348864 of
# Classic McEliece; the scheme is not compatible with Classic McEliece.
const mceliece_toy = instantiate_mceliece(6, 64, 5, "mceliece-toy")
const mceliece_light = instantiate_mceliece(10, 1024, 50, "mceliece-light")
const mceliece_standard = instantiate_mceliece(12, 3488, 64, "mceliece-standard")

const mceliece_all = [mceliece_toy, mceliece_light, mceliece_standard]

# Random invertible `k x k` matrix over GF(2).
function _random_invertible_gf2(rng::AbstractRNG, k::Int)
    while true
        S = matrix(GF(2), rand(rng, 0:1, k, k))
        is_invertible(S) && return S
    end
end

# Random binary `b x n` matrix whose rows have weight `t`.
function _random_errors(rng::AbstractRNG, b::Int, n::Int, t::Int)
    E = zeros(Int, b, n)
    for i in 1:b
        E[i, randperm(rng, n)[1:t]] .= 1
    end
    return E
end

@doc raw"""
    keygen(p::McEliece; rng::AbstractRNG = RandomDevice()) -> Tuple{FqMatrix, Tuple{BinaryGoppaCode, FqMatrix, Vector{Int}}}

Generate a key pair for the McEliece parameters `p`. A random binary Goppa
code with generator matrix $G$ is chosen together with a random invertible
$k \times k$ matrix $S$ and a random permutation matrix $P$ over
$\mathrm{GF}(2)$. The public key is the $k \times n$ matrix $SGP$, the secret
key is the triple `(C, S^{-1}, perm)` of the code, the inverse of $S$ and the
permutation. Return `(pk, sk)`.

# Examples
```jldoctest
julia> pk, sk = keygen(mceliece_toy);

julia> size(pk)
(34, 64)
```
"""
function keygen(p::McEliece; rng::AbstractRNG = RandomDevice())
    C = _goppa_code(rng, p.m, p.n, p.t)
    S = _random_invertible_gf2(rng, p.k)
    perm = randperm(rng, p.n)
    pk = (S * C.G)[:, perm]
    return pk, (C, inv(S), perm)
end

@doc raw"""
    encrypt(p::McEliece, pk::FqMatrix, msg::AbstractVector{UInt8}; rng::AbstractRNG = RandomDevice()) -> FqMatrix

Encrypt the bytes `msg` under the public key `pk`. The bits of the message are
padded with a one and zeros and split into blocks $x$ of $k$ bits, each of
which is encrypted as $x \cdot pk + e$ with a random error vector $e$ of
weight $t$. The ciphertext is the matrix over $\mathrm{GF}(2)$ with these
vectors as rows.

# Examples
```jldoctest
julia> pk, sk = keygen(mceliece_toy);

julia> ct = encrypt(mceliece_toy, pk, UInt8[1, 2, 3]);

julia> size(ct)
(1, 64)
```
"""
function encrypt(p::McEliece, pk::FqMatrix, msg::AbstractVector{UInt8}; rng::AbstractRNG = RandomDevice())
    bits = _cbc_pad(_bytes_to_bits(msg), p.k)
    b = length(bits) ÷ p.k
    X = matrix(GF(2), b, p.k, bits)
    return X * pk + matrix(GF(2), _random_errors(rng, b, p.n, p.t))
end

@doc raw"""
    decrypt(p::McEliece, sk, ct::FqMatrix) -> Union{Vector{UInt8}, Nothing}

Decrypt the ciphertext `ct` with the secret key `sk`. For each row, the
permutation is undone, the errors are removed with Patterson's algorithm
[Pat75](@cite), and the message block is obtained from the codeword with the
generator matrix of the code and the inverse of $S$. Returns `nothing` if a
row cannot be decoded or the result is not a correctly padded message.

# Examples
```jldoctest
julia> pk, sk = keygen(mceliece_toy);

julia> ct = encrypt(mceliece_toy, pk, Vector{UInt8}("oscar"));

julia> String(decrypt(mceliece_toy, sk, ct))
"oscar"
```
"""
function decrypt(p::McEliece, sk, ct::FqMatrix)
    C, S_inv, perm = sk
    @req ncols(ct) == p.n "ciphertext does not match the parameters"
    R = _gf2_int(ct[:, invperm(perm)])
    syndromes = _gf2_int(C.H * transpose(matrix(GF(2), R)))
    for i in 1:nrows(R)
        pos = _goppa_decode(C, syndromes[:, i])
        pos === nothing && return nothing
        R[i, pos] .⊻= 1
    end
    ok, X = can_solve_with_solution(C.G, matrix(GF(2), R); side = :left)
    ok || return nothing
    bits = _cbc_unpad(vec(permutedims(_gf2_int(X * S_inv))))
    return bits === nothing ? nothing : _bits_to_bytes(bits)
end
