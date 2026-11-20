# ---------------------------------------------------------------------------
# The codes used by HQC: a Reed-Solomon code over GF(2^8) concatenated with a
# duplicated first order Reed-Muller code RM(1,7). Each byte of a codeword of
# the Reed-Solomon code is encoded by RM(1,7) to 128 bits, which are repeated
# `mult` times.
# ---------------------------------------------------------------------------

# Field GF(2^8) of the Reed-Solomon code with its generator, which is a
# primitive element since the field is defined by a Conway polynomial.
_rs_field() = finite_field(2, 8, :a)

# Conversions between bytes and elements of GF(2^8).
_rs_symbol(F::FqField, a::FqFieldElem, b::UInt8) = _gf2m_element(F, a, _bytes_to_bits([b]))
_rs_byte(c::FqFieldElem) = only(_bits_to_bytes(_gf2m_bits(c)))

# Generator polynomial (x - a)(x - a^2)...(x - a^(n1 - k)) of the Reed-Solomon
# code of length `n1` and dimension `k` over the field of the generator `x`.
function _rs_generator(x::FqPolyRingElem, a::FqFieldElem, n1::Int, k::Int)
    return prod(x - a^i for i in 1:n1-k)
end

# Codeword of the Reed-Solomon code of length `n1` and dimension `k` for the
# `k` message bytes `msg`, as a vector of `n1` bytes. The message polynomial is
# multiplied by the generator polynomial.
function _rs_encode(msg::Vector{UInt8}, n1::Int)
    F, a = _rs_field()
    _, x = polynomial_ring(F, :x; cached = false)
    k = length(msg)
    m = sum(_rs_symbol(F, a, msg[i]) * x^(i - 1) for i in 1:k)
    c = m * _rs_generator(x, a, n1, k)
    return [_rs_byte(coeff(c, i)) for i in 0:n1-1]
end

# Message of `k` bytes of the received word `r` of `n1` bytes of the
# Reed-Solomon code, correcting up to (n1 - k) ÷ 2 errors. The error locator
# and evaluator are found with the extended Euclidean algorithm from the
# syndromes, the errors with Forney's formula. Returns `nothing` if there
# are too many errors.
function _rs_decode(r::Vector{UInt8}, k::Int)
    F, a = _rs_field()
    _, x = polynomial_ring(F, :x; cached = false)
    n1 = length(r)
    d = n1 - k
    c = sum(_rs_symbol(F, a, r[i]) * x^(i - 1) for i in 1:n1)
    S = sum(evaluate(c, a^j) * x^(j - 1) for j in 1:d)
    if !is_zero(S)
        r0, r1 = x^d, S
        s0, s1 = zero(x), one(x)
        while 2 * degree(r1) >= d
            q, rem = divrem(r0, r1)
            r0, r1 = r1, rem
            s0, s1 = s1, s0 - q * s1
        end
        sigma, omega = s1, r1
        dsigma = derivative(sigma)
        nerr = 0
        for i in 0:n1-1
            z = a^(-i)
            is_zero(evaluate(sigma, z)) || continue
            c += divexact(evaluate(omega, z), evaluate(dsigma, z)) * x^i
            nerr += 1
        end
        nerr == degree(sigma) || return nothing
    end
    m, rem = divrem(c, _rs_generator(x, a, n1, k))
    is_zero(rem) || return nothing
    return [_rs_byte(coeff(m, i)) for i in 0:k-1]
end

# Codeword of RM(1,7) for the byte `b`, repeated `mult` times. Bit j of a copy
# is b_0 + <(b_1, ..., b_7), j> mod 2, where j is read as a vector of 7 bits.
function _rm_encode(b::UInt8, mult::Int)
    word = [Int(b & 1) ⊻ (count_ones((b >> 1) & j) & 1) for j in 0:127]
    return repeat(word, mult)
end

# Byte whose RM(1,7) codeword repeated `mult` times is closest to the bits
# `r`. The copies are added as values +1 and -1, and the correlations with all
# codewords are computed with a fast Hadamard transform.
function _rm_decode(r::AbstractVector{Int}, mult::Int)
    f = zeros(Int, 128)
    for i in 1:mult, j in 1:128
        f[j] += 1 - 2 * r[(i - 1) * 128 + j]
    end
    h = 1
    while h < 128
        for i in 0:2h:127, j in i+1:i+h
            f[j], f[j + h] = f[j] + f[j + h], f[j] - f[j + h]
        end
        h *= 2
    end
    best = argmax(abs.(f)) - 1
    return UInt8(2 * best + (f[best + 1] < 0 ? 1 : 0))
end
