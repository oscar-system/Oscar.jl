module UOV

using Random: rand
using Keccak: shake_256
using Oscar: GF, Nemo, AbstractAlgebra

# ---------------------------------------------------------------------------
# Field setup: build fast multiplication and inverse tables using Oscar's
# finite field implementation. We use the UOV reference polynomials
#   GF(256): x^8 + x^4 + x^3 + x + 1   (0x11B)
#   GF(16) : x^4 + x + 1               (0x13)
# so that the resulting signatures match the official UOV KAT vectors.
# The tables are built once and cached per field size.
# ---------------------------------------------------------------------------

const _MUL_TAB = Dict{Int,Matrix{UInt8}}()
const _INV_TAB = Dict{Int,Vector{UInt8}}()

function _field_tables(gf::Int)
    haskey(_MUL_TAB, gf) && return _MUL_TAB[gf], _INV_TAB[gf]
    F2 = GF(2)
    R, x = AbstractAlgebra.polynomial_ring(F2, :x)
    f = (gf == 256) ? x^8 + x^4 + x^3 + x + 1 : x^4 + x + 1
    K = GF(f)
    nb = (gf == 256) ? 8 : 4
    byte2elem(b) = begin
        r = zero(K)
        for i in 0:nb-1
            if (b >> i) & 1 == 1
                r += Nemo.gen(K)^i
            end
        end
        r
    end
    elem2byte(el) = UInt8(sum([(Int(Nemo._coeff(el, i)) & 1) << i for i in 0:nb-1]))
    be = [byte2elem(b) for b in 0:gf-1]
    tab = zeros(UInt8, gf, gf)
    invtab = zeros(UInt8, gf)
    for a in 0:gf-1
        for b in 0:gf-1
            tab[a+1, b+1] = elem2byte(be[a+1] * be[b+1])
        end
    end
    for a in 1:gf-1
        invtab[a+1] = elem2byte(inv(be[a+1]))
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
# Field arithmetic (using Oscar-generated tables)
# ---------------------------------------------------------------------------

@inline gf_mul(uov::UOV, a::Integer, b::Integer) =
    uov.mul_tab[Int(a)+1, Int(b)+1]

@inline gf_inv(uov::UOV, a::Integer) = uov.inv_tab[Int(a)+1]

# multiply a vector (of m GF elements) by a scalar, elementwise
function gf_mulm(uov::UOV, vec::Vector{UInt8}, scalar::Integer)
    r = Vector{UInt8}(undef, length(vec))
    row = view(uov.mul_tab, Int(scalar)+1, :)
    @inbounds for i in eachindex(vec)
        r[i] = row[Int(vec[i])+1]
    end
    return r
end

# elementwise XOR (addition in characteristic 2), in place on dest
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
        return UInt8.(v)
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
        return UInt8.(b)
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
    b = UInt8[]
    for i in 1:d
        for j in i:d
            append!(b, pack_vec(uov, mtx[i][j]))
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
    b = UInt8[]
    for i in 1:h
        for j in 1:w
            append!(b, pack_vec(uov, mtx[i][j]))
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
            mtx[i] = UInt8.(b[p:p+uov.v-1])
        else
            mtx[i] = unpack_vec(uov, b[p:p+(uov.v÷2)-1])
        end
        p += step
    end
    return mtx
end

# Simple CTR mode (replaces AES-128-CTR)
function simple_ctr(key::Vector{UInt8}, l::Int, ctr::Int=0)::Vector{UInt8}
    result = UInt8[]
    block_ctr = ctr
    while length(result) < l
        block = zeros(UInt8, 16)
        for i in 0:15
            block[i+1] = UInt8((block_ctr >> (15 - i) * 8) & 0xFF)
        end
        aes_block = shake_256(vcat(key, block), 16)
        append!(result, aes_block)
        block_ctr += 1
    end
    return result[1:l]
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

    for j in 1:uov.m
        for i in 1:uov.v
            t = copy(m2[i][j])
            for k in i:uov.v
                vecxor!(t, gf_mulm(uov, m1[i][k], mo[j][k]))
            end
            for k in 1:uov.m
                u = gf_mulm(uov, t, mo[k][i])
                if j < k
                    vecxor!(m3[j][k], u)
                else
                    vecxor!(m3[k][j], u)
                end
            end
        end
    end

    for i in 1:uov.v
        for j in 1:uov.m
            t = copy(m2[i][j])
            for k in 1:i
                vecxor!(t, gf_mulm(uov, m1[k][i], mo[j][k]))
            end
            for k in i:uov.v
                vecxor!(t, gf_mulm(uov, m1[i][k], mo[j][k]))
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

# Apply public map to z
function pubmap(uov::UOV, z::Vector{UInt8}, tm::Vector{UInt8})
    v = uov.v
    m = uov.m

    m1 = unpack_mtri(uov, tm[1:uov.p1_sz], v)
    m2 = unpack_mrect(uov, tm[uov.p1_sz+1:uov.p1_sz+uov.p2_sz], v, m)
    m3 = unpack_mtri(uov, tm[uov.p1_sz+uov.p2_sz+1:end], m)
    x = unpack_vec(uov, z)

    y = zeros(UInt8, uov.m)
    for i in 1:v
        for j in i:v
            vecxor!(y, gf_mulm(uov, m1[i][j], gf_mul(uov, x[i], x[j])))
        end
    end
    for i in 1:v
        for j in 1:m
            vecxor!(y, gf_mulm(uov, m2[i][j], gf_mul(uov, x[i], x[v+j])))
        end
    end
    for i in 1:m
        for j in i:m
            vecxor!(y, gf_mulm(uov, m3[i][j], gf_mul(uov, x[v+i], x[v+j])))
        end
    end

    return pack_vec(uov, y)
end

# Key generation
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

# Sign a message
function sign(uov::UOV, msg::Vector{UInt8}, sk::Vector{UInt8})
    if uov.skc
        sk = expand_sk(uov, sk)
    end

    seed_sk = sk[1:uov.seed_sk_sz]
    so = sk[uov.seed_sk_sz+1:uov.seed_sk_sz+uov.so_sz]
    p1 = sk[uov.seed_sk_sz+uov.so_sz+1:uov.seed_sk_sz+uov.so_sz+uov.p1_sz]
    sks = sk[uov.seed_sk_sz+uov.so_sz+uov.p1_sz+1:end]

    m1 = unpack_mtri(uov, p1)
    ms = unpack_mrect(uov, sks, uov.v, uov.m)
    mo = unpack_rect(uov, so)

    salt = uov.rbg(uov.salt_sz)
    t = shake_256(vcat(msg, salt), uov.m_sz)

    ctr = 0
    x = nothing
    v = UInt8[]
    while x === nothing && ctr < 256
        v = unpack_vec(uov, shake_256(vcat(msg, salt, seed_sk, UInt8[ctr]), uov.v_sz))
        ctr += 1

        ll = [zeros(UInt8, uov.m) for _ in 1:uov.m]
        for i in 1:uov.m
            acc = zeros(UInt8, uov.m)
            for j in 1:uov.v
                vecxor!(acc, gf_mulm(uov, ms[j][i], v[j]))
            end
            ll[i] = acc
        end

        r = unpack_vec(uov, t)
        for i in 1:uov.v
            u = zeros(UInt8, uov.m)
            for j in i:uov.v
                vecxor!(u, gf_mulm(uov, m1[i][j], v[j]))
            end
            vecxor!(r, gf_mulm(uov, u, v[i]))
        end

        x = gauss_solve(uov, ll, r)
    end

    y = copy(v)
    for i in 1:uov.m
        for j in 1:uov.v
            y[j] ⊻= gf_mul(uov, mo[i][j], x[i])
        end
    end

    sig = vcat(pack_vec(uov, y), pack_vec(uov, x), salt)
    return sig
end

# Verify a signature
function verify(uov::UOV, sig::Vector{UInt8}, msg::Vector{UInt8}, pk::Vector{UInt8})
    if uov.pkc
        pk = expand_pk(uov, pk)
    end

    z = sig[1:uov.n_sz]
    salt = sig[uov.n_sz+1:end]

    t = shake_256(vcat(msg, salt), uov.m_sz)

    return t == pubmap(uov, z, pk)
end

# Open a signed message
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

# Instantiate UOV parameter sets
function make_uov(gf::Int, n::Int, m::Int, pkc::Bool, skc::Bool, name::String)
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
const uov_1p = make_uov(256, 112, 44, false, false, "uov-Ip-classic")
const uov_1p_pkc = make_uov(256, 112, 44, true, false, "uov-Ip-pkc")
const uov_1p_pkc_skc = make_uov(256, 112, 44, true, true, "uov-Ip-pkc-skc")
const uov_1s = make_uov(16, 160, 64, false, false, "uov-Is-classic")
const uov_1s_pkc = make_uov(16, 160, 64, true, false, "uov-Is-pkc")
const uov_1s_pkc_skc = make_uov(16, 160, 64, true, true, "uov-Is-pkc-skc")
const uov_3 = make_uov(256, 184, 72, false, false, "uov-III-classic")
const uov_3_pkc = make_uov(256, 184, 72, true, false, "uov-III-pkc")
const uov_3_pkc_skc = make_uov(256, 184, 72, true, true, "uov-III-pkc-skc")
const uov_5 = make_uov(256, 244, 96, false, false, "uov-V-classic")
const uov_5_pkc = make_uov(256, 244, 96, true, false, "uov-V-pkc")
const uov_5_pkc_skc = make_uov(256, 244, 96, true, true, "uov-V-pkc-skc")

const uov_all = [uov_1p, uov_1p_pkc, uov_1p_pkc_skc, uov_1s, uov_1s_pkc, uov_1s_pkc_skc, uov_3, uov_3_pkc, uov_3_pkc_skc, uov_5, uov_5_pkc, uov_5_pkc_skc]

# Exports
export keygen, sign, verify, open
export uov_1p, uov_1p_pkc, uov_1p_pkc_skc
export uov_1s, uov_1s_pkc, uov_1s_pkc_skc
export uov_3, uov_3_pkc, uov_3_pkc_skc
export uov_5, uov_5_pkc, uov_5_pkc_skc
export uov_all

end # module
