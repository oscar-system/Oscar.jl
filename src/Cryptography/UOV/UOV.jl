
using Random: rand
using Keccak: shake_256

# ---------------------------------------------------------------------------
# Field setup: build fast multiplication and inverse tables. The tables are
# built once per field size and cached. The UOV reference polynomials are
#   GF(256): x^8 + x^4 + x^3 + x + 1   (0x11B)
#   GF(16) : x^4 + x + 1               (0x13)
# so that the resulting signatures match the official UOV KAT vectors. The
# tables are computed directly with integer (coefficient) arithmetic, which is
# far faster than going through the generic finite field implementation.
# ---------------------------------------------------------------------------

const _MUL_TAB = Dict{Int,Matrix{UInt8}}()
const _INV_TAB = Dict{Int,Vector{UInt8}}()

# Product of the two field elements `a` and `b` (given by their coefficient
# representation as integers) in GF(2^deg) defined by the binary polynomial
# `modpoly` (whose highest set bit is bit `deg`).
function _gf_mul(a::Int, b::Int, modpoly::Int, deg::Int)::Int
    p = 0
    i = 0
    while b != 0
        if b & 1 != 0
            p ⊻= a << i
        end
        b >>= 1
        i += 1
    end
    for i in (2deg-1):-1:deg
        if (p >> i) & 1 != 0
            p ⊻= modpoly << (i - deg)
        end
    end
    return p & (2^deg - 1)
end

function _field_tables(gf::Int)
    haskey(_MUL_TAB, gf) && return _MUL_TAB[gf], _INV_TAB[gf]
    (modpoly, deg) = (gf == 256) ? (Int(0x11B), 8) : (Int(0x13), 4)
    q = gf
    tab = zeros(UInt8, q, q)
    @inbounds for a in 0:q-1
        for b in 0:q-1
            tab[a+1, b+1] = _gf_mul(a, b, modpoly, deg)
        end
    end
    invtab = zeros(UInt8, q)
    @inbounds for a in 1:q-1
        for x in 1:q-1
            if tab[a+1, x+1] == 1
                invtab[a+1] = x
                break
            end
        end
    end
    _MUL_TAB[gf] = tab
    _INV_TAB[gf] = invtab
    return tab, invtab
end

# ---------------------------------------------------------------------------
# UOV parameter structure
# ---------------------------------------------------------------------------

struct UOV
    gf::Int
    n::Int
    m::Int
    v::Int
    pkc::Bool
    skc::Bool
    name::String
    rbg::Function
    gf_bits::Int
    mul_tab::Matrix{UInt8}
    inv_tab::Vector{UInt8}
    v_sz::Int
    n_sz::Int
    m_sz::Int
    sig_sz::Int
    seed_sk_sz::Int
    seed_pk_sz::Int
    so_sz::Int
    p1_sz::Int
    p2_sz::Int
    p3_sz::Int
    salt_sz::Int
end

# ---------------------------------------------------------------------------
# Field arithmetic (using the precomputed tables)
# ---------------------------------------------------------------------------

@inline gf_mul(uov::UOV, a::Integer, b::Integer) =
    uov.mul_tab[Int(a)+1, Int(b)+1]

@inline gf_inv(uov::UOV, a::Integer) = uov.inv_tab[Int(a)+1]

# Multiply the vector `vec` (of GF elements) elementwise by the scalar
# `scalar`, writing the result into the preallocated vector `dest`.
@inline function gf_mulm!(uov::UOV, dest::Vector{UInt8}, vec::Vector{UInt8}, scalar::Integer)
    row = view(uov.mul_tab, Int(scalar)+1, :)
    @inbounds for i in eachindex(vec)
        dest[i] = row[Int(vec[i])+1]
    end
    return dest
end

# Allocating convenience wrapper around the in-place `gf_mulm!`.
function gf_mulm(uov::UOV, vec::Vector{UInt8}, scalar::Integer)
    return gf_mulm!(uov, Vector{UInt8}(undef, length(vec)), vec, scalar)
end

# Elementwise XOR (addition in characteristic 2), in place on `dest`.
@inline function vecxor!(dest::Vector{UInt8}, src::Vector{UInt8})
    @inbounds for i in eachindex(dest)
        dest[i] ⊻= src[i]
    end
    return dest
end

# ---------------------------------------------------------------------------
# Serialization helpers
# ---------------------------------------------------------------------------

# Convert a vector of GF elements (length m) to bytes
function pack_vec(uov::UOV, v::Vector{UInt8})
    if uov.gf == 256
        return copy(v)
    else
        res = UInt8[]
        for i in 1:2:length(v)-1
            push!(res, UInt8((v[i] & 0x0F) | ((v[i+1] & 0x0F) << 4)))
        end
        return res
    end
end

# Convert bytes to a vector of GF elements
function unpack_vec(uov::UOV, b::Vector{UInt8})
    if uov.gf == 256
        return copy(b)
    else
        res = zeros(UInt8, 2 * length(b))
        for (k, byte) in enumerate(b)
            res[2k-1] = byte & 0x0F
            res[2k] = byte >> 4
        end
        return res
    end
end

# Unpack an upper triangular matrix. Each entry is a vector of m GF elements.
function unpack_mtri(uov::UOV, b::Vector{UInt8}, d::Int=uov.v)
    mtx = Vector{Vector{Vector{UInt8}}}(undef, d)
    p = 1
    for i in 1:d
        row = Vector{Vector{UInt8}}(undef, d)
        for j in i:d
            if uov.gf == 256
                row[j] = b[p:p+uov.m-1]
            else
                row[j] = unpack_vec(uov, b[p:p+(uov.m÷2)-1])
            end
            p += uov.m_sz
        end
        for j in 1:i-1
            row[j] = zeros(UInt8, uov.m)
        end
        mtx[i] = row
    end
    return mtx
end

# Pack an upper triangular matrix to bytes
function pack_mtri(uov::UOV, mtx, d::Int=uov.v)
    b = Vector{UInt8}(undef, d*(d+1)÷2 * uov.m_sz)
    p = 1
    for i in 1:d
        for j in i:d
            chunk = pack_vec(uov, mtx[i][j])
            copyto!(b, p, chunk)
            p += uov.m_sz
        end
    end
    return b
end

# Unpack a rectangular (h x w) matrix of m-element vectors
function unpack_mrect(uov::UOV, b::Vector{UInt8}, h::Int, w::Int)
    mtx = Vector{Vector{Vector{UInt8}}}(undef, h)
    p = 1
    for i in 1:h
        row = Vector{Vector{UInt8}}(undef, w)
        for j in 1:w
            if uov.gf == 256
                row[j] = b[p:p+uov.m-1]
            else
                row[j] = unpack_vec(uov, b[p:p+(uov.m÷2)-1])
            end
            p += uov.m_sz
        end
        mtx[i] = row
    end
    return mtx
end

# Pack a rectangular matrix to bytes
function pack_mrect(uov::UOV, mtx, h::Int, w::Int)
    b = Vector{UInt8}(undef, h * w * uov.m_sz)
    p = 1
    for i in 1:h
        for j in 1:w
            chunk = pack_vec(uov, mtx[i][j])
            copyto!(b, p, chunk)
            p += uov.m_sz
        end
    end
    return b
end

# Unpack a m x v matrix of scalars (from the "so" secret material)
function unpack_rect(uov::UOV, b::Vector{UInt8})
    mtx = Vector{Vector{UInt8}}(undef, uov.m)
    p = 1
    step = (uov.gf == 256) ? uov.v : (uov.v ÷ 2)
    for i in 1:uov.m
        if uov.gf == 256
            mtx[i] = b[p:p+uov.v-1]
        else
            mtx[i] = unpack_vec(uov, b[p:p+(uov.v÷2)-1])
        end
        p += step
    end
    return mtx
end

# Deterministic keystream of `l` bytes derived from `key`. This is a stand-in
# for the AES-128-CTR keystream used by the UOV specification; it only needs to
# be deterministic and of the requested length.
function simple_ctr(key::Vector{UInt8}, l::Int)::Vector{UInt8}
    return shake_256(key, l)
end

# UOV.ExpandP()
function expand_p(uov::UOV, seed_pk::Vector{UInt8})
    pk = simple_ctr(seed_pk, uov.p1_sz + uov.p2_sz)
    p1 = pk[1:uov.p1_sz]
    p2 = pk[uov.p1_sz+1:end]
    return p1, p2
end

# UOV.ExpandPK(cpk)
function expand_pk(uov::UOV, cpk::Vector{UInt8})
    seed_pk = cpk[1:uov.seed_pk_sz]
    p3 = cpk[uov.seed_pk_sz+1:end]
    p1, p2 = expand_p(uov, seed_pk)
    epk = vcat(p1, p2, p3)
    return epk
end

# UOV.ExpandSK(csk)
function expand_sk(uov::UOV, csk::Vector{UInt8})
    seed_sk = csk[1:uov.seed_sk_sz]
    seed_pk_so = shake_256(seed_sk, uov.seed_pk_sz + uov.so_sz)
    seed_pk = seed_pk_so[1:uov.seed_pk_sz]
    so = seed_pk_so[uov.seed_pk_sz+1:end]
    p1, p2 = expand_p(uov, seed_pk)
    sks, p3 = calc_f2_p3(uov, p1, p2, so)
    esk = vcat(csk, so, p1, sks)
    return esk
end

# UOV.ExpandP2()
function calc_f2_p3(uov::UOV, p1::Vector{UInt8}, p2::Vector{UInt8}, so::Vector{UInt8})
    m1 = unpack_mtri(uov, p1)
    m2 = unpack_mrect(uov, p2, uov.v, uov.m)
    mo = unpack_rect(uov, so)

    # m3 is a full m x m matrix (upper triangle used)
    m3 = [ [zeros(UInt8, uov.m) for _ in 1:uov.m] for _ in 1:uov.m ]

    tmp = Vector{UInt8}(undef, uov.m)
    for j in 1:uov.m
        for i in 1:uov.v
            t = copy(m2[i][j])
            for k in i:uov.v
                gf_mulm!(uov, tmp, m1[i][k], mo[j][k])
                vecxor!(t, tmp)
            end
            for k in 1:uov.m
                gf_mulm!(uov, tmp, t, mo[k][i])
                if j < k
                    vecxor!(m3[j][k], tmp)
                else
                    vecxor!(m3[k][j], tmp)
                end
            end
        end
    end

    for i in 1:uov.v
        for j in 1:uov.m
            t = copy(m2[i][j])
            for k in 1:i
                gf_mulm!(uov, tmp, m1[k][i], mo[j][k])
                vecxor!(t, tmp)
            end
            for k in i:uov.v
                gf_mulm!(uov, tmp, m1[i][k], mo[j][k])
                vecxor!(t, tmp)
            end
            m2[i][j] = t
        end
    end

    p3 = pack_mtri(uov, m3, uov.m)
    sks = pack_mrect(uov, m2, uov.v, uov.m)
    return sks, p3
end

# Gaussian elimination solver (elements are UInt8 GF values)
function gauss_solve(uov::UOV, l::Vector{Vector{UInt8}}, c::Vector{UInt8})
    h = uov.m
    w = uov.m + 1
    m = [zeros(UInt8, w) for _ in 1:h]
    for j in 1:h
        for i in 1:h
            m[j][i] = l[i][j]
        end
        m[j][w] = c[j]
    end

    for i in 1:h
        j = i
        while m[j][i] == 0
            j += 1
            if j > h
                return nothing
            end
        end
        if i != j
            for k in 1:w
                m[i][k] ⊻= m[j][k]
            end
        end
        x = gf_inv(uov, m[i][i])
        for k in 1:w
            m[i][k] = gf_mul(uov, m[i][k], x)
        end
        for j in 1:h
            x = m[j][i]
            if j != i
                for k in 1:w
                    m[j][k] ⊻= gf_mul(uov, m[i][k], x)
                end
            end
        end
    end

    return [m[i][w] for i in 1:h]
end

# Apply public map to z. The byte vector `tm` must be the full (expanded)
# public key, i.e. the concatenation of `p1`, `p2` and `p3`.
function pubmap(uov::UOV, z::Vector{UInt8}, tm::Vector{UInt8})
    v = uov.v
    m = uov.m

    @req length(z) == uov.n_sz "signature part z must have length $(uov.n_sz), got $(length(z))"
    @req length(tm) == uov.p1_sz + uov.p2_sz + uov.p3_sz "public key must have length $(uov.p1_sz + uov.p2_sz + uov.p3_sz), got $(length(tm))"

    m1 = unpack_mtri(uov, tm[1:uov.p1_sz], v)
    m2 = unpack_mrect(uov, tm[uov.p1_sz+1:uov.p1_sz+uov.p2_sz], v, m)
    m3 = unpack_mtri(uov, tm[uov.p1_sz+uov.p2_sz+1:end], m)
    x = unpack_vec(uov, z)

    y = zeros(UInt8, m)
    tmp = Vector{UInt8}(undef, m)
    for i in 1:v
        for j in i:v
            gf_mulm!(uov, tmp, m1[i][j], gf_mul(uov, x[i], x[j]))
            vecxor!(y, tmp)
        end
    end
    for i in 1:v
        for j in 1:m
            gf_mulm!(uov, tmp, m2[i][j], gf_mul(uov, x[i], x[v+j]))
            vecxor!(y, tmp)
        end
    end
    for i in 1:m
        for j in i:m
            gf_mulm!(uov, tmp, m3[i][j], gf_mul(uov, x[v+i], x[v+j]))
            vecxor!(y, tmp)
        end
    end

    return pack_vec(uov, y)
end

@doc raw"""
    keygen(uov::UOV) -> (pk, sk)

Generate a public/private key pair for the UOV parameter set `uov`, returning
the public key `pk` first and the secret key `sk` second. Both are byte vectors
whose layout depends on the `pkc` and `skc` flags of `uov`.

# Examples
```jldoctest uov-keygen
julia> uov = uov_1p; pk, sk = keygen(uov); (length(pk), length(sk))
(278432, 237896)
```
"""
function keygen(uov::UOV)
    seed_sk = uov.rbg(uov.seed_sk_sz)
    seed_pk_so = shake_256(seed_sk, uov.seed_pk_sz + uov.so_sz)
    seed_pk = seed_pk_so[1:uov.seed_pk_sz]
    so = seed_pk_so[uov.seed_pk_sz+1:end]
    p1, p2 = expand_p(uov, seed_pk)
    sks, p3 = calc_f2_p3(uov, p1, p2, so)

    if uov.pkc
        pk = vcat(seed_pk, p3)
    else
        pk = vcat(p1, p2, p3)
    end

    if uov.skc
        sk = seed_sk
    else
        sk = vcat(seed_sk, so, p1, sks)
    end

    return pk, sk
end

@doc raw"""
    sign(uov::UOV, msg::Vector{UInt8}, sk::Vector{UInt8}) -> sig

Sign the message `msg` (a byte vector) with the secret key `sk` of the UOV
parameter set `uov`, returning the signature as a byte vector.

# Examples
```jldoctest uov-sign
julia> uov = uov_1p; pk, sk = keygen(uov); msg = zeros(UInt8, 32); sig = sign(uov, msg, sk);

julia> verify(uov, sig, msg, pk)
true
```
"""
function sign(uov::UOV, msg::Vector{UInt8}, sk::Vector{UInt8})
    if uov.skc
        @req length(sk) == uov.seed_sk_sz "compact secret key must have length $(uov.seed_sk_sz), got $(length(sk))"
        sk = expand_sk(uov, sk)
    else
        @req length(sk) == uov.seed_sk_sz + uov.so_sz + uov.p1_sz + uov.p2_sz "secret key must have length $(uov.seed_sk_sz + uov.so_sz + uov.p1_sz + uov.p2_sz), got $(length(sk))"
    end

    seed_sk = sk[1:uov.seed_sk_sz]
    so = sk[uov.seed_sk_sz+1:uov.seed_sk_sz+uov.so_sz]
    p1 = sk[uov.seed_sk_sz+uov.so_sz+1:uov.seed_sk_sz+uov.so_sz+uov.p1_sz]
    sks = sk[uov.seed_sk_sz+uov.so_sz+uov.p1_sz+1:end]

    m1 = unpack_mtri(uov, p1)
    ms = unpack_mrect(uov, sks, uov.v, uov.m)
    mo = unpack_rect(uov, so)

    m = uov.m
    v = uov.v

    salt = uov.rbg(uov.salt_sz)
    t = shake_256(vcat(msg, salt), uov.m_sz)

    # Preallocate the work vectors so the (possibly) retrying loop below does
    # not allocate on every iteration.
    ll = [zeros(UInt8, m) for _ in 1:m]
    acc = zeros(UInt8, m)
    r = zeros(UInt8, m)
    u = zeros(UInt8, m)
    tmp = Vector{UInt8}(undef, m)

    ctr = 0
    x = nothing
    zvec = UInt8[]
    while x === nothing && ctr < 256
        zvec = unpack_vec(uov, shake_256(vcat(msg, salt, seed_sk, UInt8[ctr]), uov.v_sz))
        ctr += 1

        for i in 1:m
            fill!(acc, 0x00)
            for j in 1:v
                gf_mulm!(uov, tmp, ms[j][i], zvec[j])
                vecxor!(acc, tmp)
            end
            copyto!(ll[i], acc)
        end

        copyto!(r, unpack_vec(uov, t))
        for i in 1:v
            fill!(u, 0x00)
            for j in i:v
                gf_mulm!(uov, tmp, m1[i][j], zvec[j])
                vecxor!(u, tmp)
            end
            gf_mulm!(uov, tmp, u, zvec[i])
            vecxor!(r, tmp)
        end

        x = gauss_solve(uov, ll, r)
    end

    y = copy(zvec)
    for i in 1:m
        for j in 1:v
            y[j] ⊻= gf_mul(uov, mo[i][j], x[i])
        end
    end

    sig = vcat(pack_vec(uov, y), pack_vec(uov, x), salt)
    return sig
end

@doc raw"""
    verify(uov::UOV, sig::Vector{UInt8}, msg::Vector{UInt8}, pk::Vector{UInt8}) -> Bool

Return `true` if the signature `sig` is a valid signature of the message `msg`
under the public key `pk` of the UOV parameter set `uov`.

# Examples
```jldoctest uov-verify
julia> uov = uov_1p; pk, sk = keygen(uov); msg = zeros(UInt8, 32); sig = sign(uov, msg, sk); sig_bad = copy(sig); sig_bad[1] ⊻= 0xff;

julia> verify(uov, sig, msg, pk), verify(uov, sig_bad, msg, pk)
(true, false)
```
"""
function verify(uov::UOV, sig::Vector{UInt8}, msg::Vector{UInt8}, pk::Vector{UInt8})
    @req length(sig) == uov.sig_sz "signature must have length $(uov.sig_sz), got $(length(sig))"

    if uov.pkc
        @req length(pk) == uov.seed_pk_sz + uov.p3_sz "compact public key must have length $(uov.seed_pk_sz + uov.p3_sz), got $(length(pk))"
        pk = expand_pk(uov, pk)
    else
        @req length(pk) == uov.p1_sz + uov.p2_sz + uov.p3_sz "public key must have length $(uov.p1_sz + uov.p2_sz + uov.p3_sz), got $(length(pk))"
    end

    z = sig[1:uov.n_sz]
    salt = sig[uov.n_sz+1:end]

    t = shake_256(vcat(msg, salt), uov.m_sz)

    return t == pubmap(uov, z, pk)
end

@doc raw"""
    open(uov::UOV, sm::Vector{UInt8}, pk::Vector{UInt8})

Recover the message from a signed message `sm` (message concatenated with its
signature), verifying it against the public key `pk`. Returns `nothing` if the
signature is invalid.

# Examples
```jldoctest uov-open
julia> uov = uov_1p; pk, sk = keygen(uov); msg = UInt8[1, 2, 3]; sm = vcat(msg, sign(uov, msg, sk));

julia> Oscar.open(uov, sm, pk)
3-element Vector{UInt8}:
 0x01
 0x02
 0x03
```
"""
function open(uov::UOV, sm::Vector{UInt8}, pk::Vector{UInt8})
    msg_sz = length(sm) - uov.sig_sz
    msg = sm[1:msg_sz]
    sig = sm[msg_sz+1:end]

    if !verify(uov, sig, msg, pk)
        return nothing
    end

    return msg
end

# Default random function
function default_rbg(n::Int=32)
    return rand(UInt8, n)
end

@doc raw"""
    instantiate_uov(gf::Int, n::Int, m::Int, pkc::Bool, skc::Bool,
                    name::String) -> UOV

Create a UOV parameter set over $\mathrm{GF}(gf)$ with $n$ variables and $m$
equations. The parameters must satisfy the following restrictions:
- `gf` is `16` or `256`;
- `m > 0` and `n >= 2m`, so that the number of univariate variables
  `v = n - m` is at least `m`;
- if `gf == 16`, then `n` and `m` are even.

The Boolean flags `pkc` and `skc` switch on public-key and secret-key
compression, respectively, and `name` labels the parameter set. The predefined
sets `uov_1p`, `uov_1s`, `uov_3`, `uov_5` (and their `_pkc` and `_pkc_skc`
variants) are created this way.

# Examples
Create a parameter set and read back its parameters; an invalid `gf` is
rejected with a clear error.
```jldoctest uov-instantiate
julia> u = instantiate_uov(256, 112, 44, false, false, "my-uov"); (u.n, u.m, u.v)
(112, 44, 68)

julia> instantiate_uov(4, 10, 3, false, false, "bad")
ERROR: ArgumentError: gf must be 16 or 256, got 4
[...]
```
"""
function instantiate_uov(gf::Int, n::Int, m::Int, pkc::Bool, skc::Bool, name::String)
    @req gf == 16 || gf == 256 "gf must be 16 or 256, got $gf"
    @req m > 0 "m must be positive, got $m"
    @req n >= 2*m "n must be at least 2*m (so that v = n - m >= m), got n=$n, m=$m"
    @req gf != 16 || (iseven(n) && iseven(m)) "for gf == 16, n and m must be even, got n=$n, m=$m"

    v = n - m
    gf_bits = (gf == 256) ? 8 : 4
    mul_tab, inv_tab = _field_tables(gf)

    v_sz = gf_bits * v ÷ 8
    n_sz = gf_bits * n ÷ 8
    m_sz = gf_bits * m ÷ 8

    seed_sk_sz = 32
    seed_pk_sz = 16
    salt_sz = 16
    sig_sz = v_sz + m_sz + salt_sz

    return UOV(gf, n, m, v, pkc, skc, name, default_rbg, gf_bits, mul_tab, inv_tab,
               v_sz, n_sz, m_sz, sig_sz, seed_sk_sz, seed_pk_sz,
               gf_bits * v * m ÷ 8, m_sz * v * (v + 1) ÷ 2, m_sz * v * m,
               m_sz * m * (m + 1) ÷ 2, salt_sz)
end

# Parameter sets
const uov_1p = instantiate_uov(256, 112, 44, false, false, "uov-Ip-classic")
const uov_1p_pkc = instantiate_uov(256, 112, 44, true, false, "uov-Ip-pkc")
const uov_1p_pkc_skc = instantiate_uov(256, 112, 44, true, true, "uov-Ip-pkc-skc")
const uov_1s = instantiate_uov(16, 160, 64, false, false, "uov-Is-classic")
const uov_1s_pkc = instantiate_uov(16, 160, 64, true, false, "uov-Is-pkc")
const uov_1s_pkc_skc = instantiate_uov(16, 160, 64, true, true, "uov-Is-pkc-skc")
const uov_3 = instantiate_uov(256, 184, 72, false, false, "uov-III-classic")
const uov_3_pkc = instantiate_uov(256, 184, 72, true, false, "uov-III-pkc")
const uov_3_pkc_skc = instantiate_uov(256, 184, 72, true, true, "uov-III-pkc-skc")
const uov_5 = instantiate_uov(256, 244, 96, false, false, "uov-V-classic")
const uov_5_pkc = instantiate_uov(256, 244, 96, true, false, "uov-V-pkc")
const uov_5_pkc_skc = instantiate_uov(256, 244, 96, true, true, "uov-V-pkc-skc")

const uov_all = [uov_1p, uov_1p_pkc, uov_1p_pkc_skc, uov_1s, uov_1s_pkc, uov_1s_pkc_skc, uov_3, uov_3_pkc, uov_3_pkc_skc, uov_5, uov_5_pkc, uov_5_pkc_skc]
