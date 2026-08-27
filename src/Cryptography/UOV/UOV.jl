module UOV

using Random: rand
using Keccak: shake_256

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
    gf_mul::Function
    gf_mulm::Function
    v_sz::Int
    n_sz::Int
    m_sz::Int
    seed_sk_sz::Int
    seed_pk_sz::Int
    so_sz::Int
    p1_sz::Int
    p2_sz::Int
    p3_sz::Int
    salt_sz::Int
    mm::BigInt
end

# GF(16) multiplication
function gf16_mul(a::Int, b::Int)::Int
    r = a & (-(b & 1))
    for i in 1:3
        t = a & 8
        a = ((a ⊻ t) << 1) ⊻ (t >> 2) ⊻ (t >> 3)
        r ⊻= a & (-(b >> i) & 1)
    end
    return r
end

# GF(16) vector * scalar multiply (BigInt version)
function gf16_mulm(v::BigInt, a::Int, mm::BigInt)::BigInt
    r = v & (-(a & 1))
    for i in 1:3
        t = v & mm
        v = ((v ⊻ t) << 1) ⊻ (t >> 3) ⊻ (t >> 2)
        if (a >> i) & 1 == 1
            r ⊻= v
        end
    end
    return r
end

# GF(16) vector * scalar multiply (Int version)
function gf16_mulm(v::Int, a::Int, mm::Int)::Int
    r = v & (-(a & 1))
    for i in 1:3
        t = v & mm
        v = ((v ⊻ t) << 1) ⊻ (t >> 3) ⊻ (t >> 2)
        if (a >> i) & 1 == 1
            r ⊻= v
        end
    end
    return r
end

# GF(256) multiplication
function gf256_mul(a::Int, b::Int)::Int
    r = a & (-(b & 1))
    for i in 1:7
        a = (a << 1) ⊻ ((-(a >> 7)) & 0x11B)
        r ⊻= a & (-(b >> i) & 1)
    end
    return r
end

# GF(256) vector * scalar multiply (BigInt version)
function gf256_mulm(v::BigInt, a::Int, mm::BigInt)::BigInt
    r = v & (-(a & 1))
    for i in 1:7
        t = v & mm
        v = ((v ⊻ t) << 1) ⊻ (t >> 7) ⊻ (t >> 6) ⊻ (t >> 4) ⊻ (t >> 3)
        if (a >> i) & 1 == 1
            r ⊻= v
        end
    end
    return r
end

# GF(256) vector * scalar multiply (Int version)
function gf256_mulm(v::Int, a::Int, mm::Int)::Int
    r = v & (-(a & 1))
    for i in 1:7
        t = v & mm
        v = ((v ⊻ t) << 1) ⊻ (t >> 7) ⊻ (t >> 6) ⊻ (t >> 4) ⊻ (t >> 3)
        if (a >> i) & 1 == 1
            r ⊻= v
        end
    end
    return r
end

# GF inverse (using extended Euclidean algorithm)
function gf_inv(a::Int, gf::Int, gf_bits::Int, gf_mul::Function)::Int
    r = a
    for _ in 2:gf_bits-1
        a = gf_mul(a, a)
        r = gf_mul(r, a)
    end
    r = gf_mul(r, r)
    return r
end

# Pack a vector of GF elements into bytes
function gf_pack(v::Vector{Int}, gf::Int)::Vector{UInt8}
    if gf == 256
        return UInt8.(v .& 0xFF)
    elseif gf == 16
        result = UInt8[]
        for i in 1:2:length(v)-1
            push!(result, UInt8((v[i] & 0xF) + ((v[i + 1] & 0xF) << 4)))
        end
        return result
    end
end

# Unpack bytes into a vector of GF elements
function gf_unpack(b::Vector{UInt8}, gf::Int)::Vector{Int}
    if gf == 256
        return Int.(b)
    elseif gf == 16
        v = Int[]
        for x in b
            push!(v, Int(x & 0xF))
            push!(v, Int(x >> 4))
        end
        return v
    end
end

# Unpack an upper triangular matrix from bytes
function unpack_mtri(b::Vector{UInt8}, d::Int, m_sz::Int)
    m = Vector{BigInt}[]
    p = 1
    for i in 1:d
        row = BigInt[]
        for j in i:d
            t = BigInt(0)
            for k in 0:m_sz-1
                t = t * 256 + BigInt(b[p + k])
            end
            push!(row, t)
            p += m_sz
        end
        push!(m, row)
    end
    return m
end

# Pack an upper triangular matrix to bytes
function pack_mtri(m::Vector{Vector{BigInt}}, d::Int, m_sz::Int)::Vector{UInt8}
    b = UInt8[]
    for i in 1:d
        for j in i:d
            t = m[i][j - i + 1]
            for k in m_sz:-1:1
                push!(b, UInt8(t & 0xFF))
                t = t >>> 8
            end
        end
    end
    return b
end

# Unpack a rectangular matrix from bytes
function unpack_mrect(b::Vector{UInt8}, h::Int, w::Int, m_sz::Int)
    m = Vector{BigInt}[]
    p = 1
    for i in 1:h
        row = BigInt[]
        for j in 1:w
            t = BigInt(0)
            for k in 0:m_sz-1
                t = t * 256 + BigInt(b[p + k])
            end
            push!(row, t)
            p += m_sz
        end
        push!(m, row)
    end
    return m
end

# Pack a rectangular matrix to bytes
function pack_mrect(m::Vector{Vector{BigInt}}, h::Int, w::Int, m_sz::Int)::Vector{UInt8}
    b = UInt8[]
    for i in 1:h
        for j in 1:w
            t = m[i][j]
            for k in m_sz:-1:1
                push!(b, UInt8(t & 0xFF))
                t = t >>> 8
            end
        end
    end
    return b
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
    m1 = unpack_mtri(p1, uov.v, uov.m_sz)
    m2 = unpack_mrect(p2, uov.v, uov.m, uov.m_sz)
    mo = [gf_unpack(so[i:i+uov.v-1], uov.gf) for i in 1:uov.v:uov.so_sz]
    
    m3 = [BigInt[BigInt(0) for _ in 1:uov.m] for _ in 1:uov.m]
    
    for j in 1:uov.m
        for i in 1:uov.v
            t = m2[i][j]
            for k in i:uov.v
                t = t ⊻ uov.gf_mulm(m1[i][k - i + 1], mo[j][k], uov.mm)
            end
            for k in 1:uov.m
                u = uov.gf_mulm(t, mo[k][i], uov.mm)
                if j < k
                    m3[j][k] = m3[j][k] ⊻ u
                else
                    m3[k][j] = m3[k][j] ⊻ u
                end
            end
        end
    end
    
    for i in 1:uov.v
        for j in 1:uov.m
            t = m2[i][j]
            for k in 1:i
                t = t ⊻ uov.gf_mulm(m1[k][i - k + 1], mo[j][k], uov.mm)
            end
            for k in i:uov.v
                t = t ⊻ uov.gf_mulm(m1[i][k - i + 1], mo[j][k], uov.mm)
            end
            m2[i][j] = t
        end
    end
    
    p3 = pack_mtri(m3, uov.m, uov.m_sz)
    sks = pack_mrect(m2, uov.v, uov.m, uov.m_sz)
    
    return sks, p3
end

# Gaussian elimination solver
function gauss_solve(uov::UOV, l::Vector{BigInt}, c::Vector{Int})
    h = uov.m
    w = uov.m + 1
    
    L = Vector{Int}[]
    for i in 1:uov.m
        bytes = int_to_bytes(l[i], uov.m_sz)
        li = gf_unpack(bytes, uov.gf)
        push!(L, li)
    end
    
    @assert length(L) == uov.m "L should have $uov.m rows, has $(length(L))"
    @assert length(L[1]) == uov.m "L[1] should have $uov.m columns, has $(length(L[1]))"
    
    m = [zeros(Int, w) for _ in 1:uov.m]
    for j in 1:uov.m
        for i in 1:uov.m
            m[j][i] = L[i][j]
        end
        m[j][w] = c[j]
    end
    
    for i in 1:uov.m
        j = i
        while m[j][i] == 0
            j += 1
            if j > uov.m
                return nothing
            end
        end
        if i != j
            for k in 1:w
                m[i][k] ⊻= m[j][k]
            end
        end
        x = gf_inv(m[i][i], uov.gf, uov.gf_bits, uov.gf_mul)
        for k in 1:w
            m[i][k] = uov.gf_mul(m[i][k], x)
        end
        for j in 1:uov.m
            x = m[j][i]
            if j != i
                for k in 1:w
                    m[j][k] ⊻= uov.gf_mul(m[i][k], x)
                end
            end
        end
    end
    
    return [m[i][w] for i in 1:uov.m]
end

# Apply public map to z
function pubmap(uov::UOV, z::Vector{UInt8}, tm::Vector{UInt8})
    v = uov.v
    m = uov.m
    
    m1 = unpack_mtri(tm[1:uov.p1_sz], v, uov.m_sz)
    m2 = unpack_mrect(tm[uov.p1_sz+1:uov.p1_sz+uov.p2_sz], v, m, uov.m_sz)
    m3 = unpack_mtri(tm[uov.p1_sz+uov.p2_sz+1:end], m, uov.m_sz)
    x = gf_unpack(z, uov.gf)
    
    y = BigInt(0)
    # P1
    for i in 1:v
        for j in i:v
            y ⊻= uov.gf_mulm(m1[i][j - i + 1], uov.gf_mul(Int(x[i]), Int(x[j])), uov.mm)
        end
    end
    
    for i in 1:v
        for j in 1:m
            y ⊻= uov.gf_mulm(m2[i][j], uov.gf_mul(Int(x[i]), Int(x[v + j])), uov.mm)
        end
    end
    
    for i in 1:m
        for j in i:m
            y ⊻= uov.gf_mulm(m3[i][j - i + 1], uov.gf_mul(Int(x[v + i]), Int(x[v + j])), uov.mm)
        end
    end
    
    return int_to_bytes(y, uov.m_sz)
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
        pk = vcat(p1, sks, p3)
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
    
    m1 = unpack_mtri(p1, uov.v, uov.m_sz)
    ms = unpack_mrect(sks, uov.v, uov.m, uov.m_sz)
    
    salt = uov.rbg(uov.salt_sz)
    t = shake_256(vcat(msg, salt), uov.m_sz)
    
    ctr = 0
    x = nothing
    v = Int[]
    while x === nothing && ctr < 256
        v = gf_unpack(shake_256(vcat(msg, salt, seed_sk, UInt8[ctr]), uov.v_sz), uov.gf)
        ctr += 1
        
        ll = [BigInt(0) for _ in 1:uov.m]
        
        for i in 1:uov.m
            for j in 1:uov.v
                ll[i] = ll[i] ⊻ uov.gf_mulm(ms[j][i], v[j], uov.mm)
            end
        end
        
        r = BigInt(0)
        for (idx, byte_val) in enumerate(t)
            r = r * 256 + BigInt(byte_val)
        end
        for i in 1:uov.v
            u = BigInt(0)
            for j in i:uov.v
                u = u ⊻ uov.gf_mulm(m1[i][j - i + 1], v[j], uov.mm)
            end
            r = r ⊻ uov.gf_mulm(u, v[i], uov.mm)
        end
        r = gf_unpack(int_to_bytes(r, uov.m_sz), uov.gf)
        
        x = gauss_solve(uov, ll, r)
        
    end
    
    y = vcat(v)
    for i in 1:uov.m
        for j in 1:uov.v
            y[j] ⊻= uov.gf_mul(Int(ms[j][i] & 0xFF), x[i])
        end
    end
    
    sig = vcat(gf_pack(y, uov.gf), gf_pack(x, uov.gf), salt)
    
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

# Int to bytes conversion
function int_to_bytes(val::Int, n::Int)::Vector{UInt8}
    result = UInt8[]
    for i in 1:n
        push!(result, UInt8((val >>> (8 * (n - i))) & 0xFF))
    end
    return result
end

function int_to_bytes(val::BigInt, n::Int)::Vector{UInt8}
    result = UInt8[]
    for i in 1:n
        push!(result, UInt8((val >>> (8 * (n - i))) & 0xFF))
    end
    return result
end

# Instantiate UOV parameter sets
function make_uov(gf::Int, n::Int, m::Int, pkc::Bool, skc::Bool, name::String)
    v = n - m
    
    if pkc
        kc = "pkc"
    else
        kc = "classic"
    end
    if skc
        kc *= "-skc"
    end
    
    katname = "OV($gf,$n,$m)-$kc"
    
    if gf == 256
        gf_bits = 8
        gf_mul = gf256_mul
        gf_mulm = gf256_mulm
    elseif gf == 16
        gf_bits = 4
        gf_mul = gf16_mul
        gf_mulm = gf16_mulm
    else
        throw(ArgumentError("Invalid gf: $gf. Must be 256 or 16."))
    end
    
    v_sz = gf_bits * v ÷ 8
    n_sz = gf_bits * n ÷ 8
    m_sz = gf_bits * m ÷ 8
    
    seed_sk_sz = 32
    seed_pk_sz = 16
    salt_sz = 16
    
    mm = BigInt(0)
    for i in 1:m
        mm = (mm << 1) ⊻ (i <= m ? 0 : 0)
    end
    if gf == 256
        mm = mm ⊻ 0x11B
    elseif gf == 16
        mm = mm ⊻ 0x13
    end
    
    return UOV(gf, n, m, v, pkc, skc, name, default_rbg, gf_bits, gf_mul, gf_mulm, v_sz, n_sz, m_sz, seed_sk_sz, seed_pk_sz, gf_bits * v * m ÷ 8, m_sz * v * (v + 1) ÷ 2, m_sz * v * m, m_sz * m * (m + 1) ÷ 2, salt_sz, mm)
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
