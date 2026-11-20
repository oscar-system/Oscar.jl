# ---------------------------------------------------------------------------
# Binary Goppa codes, the codes used by the encryption schemes of McEliece and
# Niederreiter. A binary Goppa code of length n is given by a support
# L = (a_1, ..., a_n) of distinct elements of GF(2^m) and an irreducible
# polynomial g of degree t over GF(2^m). It consists of all c in GF(2)^n with
#
#   sum_i c_i / (x - a_i) = 0 mod g,
#
# has dimension at least n - mt and corrects t errors with Patterson's
# algorithm [Pat75](@cite).
# ---------------------------------------------------------------------------

struct BinaryGoppaCode
    m::Int
    t::Int
    support::Vector{FqFieldElem}
    g::FqPolyRingElem
    H::FqMatrix  # binary parity check matrix with m*t rows and n columns
    G::FqMatrix  # binary generator matrix with n - m*t rows and n columns
end

# Coordinates of an element of GF(2^m) over GF(2), as integers 0 and 1.
_gf2m_bits(a::FqFieldElem) = [is_one(c) ? 1 : 0 for c in absolute_coordinates(a)]

# Element of GF(2^m) with the coordinates `bits` with respect to the powers of `a`.
_gf2m_element(F::FqField, a::FqFieldElem, bits::AbstractVector{Int}) =
    sum(bits[j] * a^(j - 1) for j in eachindex(bits); init = zero(F))

# Random monic irreducible polynomial of degree `t` over the field of `x`.
function _random_irreducible(rng::AbstractRNG, x::FqPolyRingElem, t::Int)
    F = base_ring(parent(x))
    while true
        g = x^t + sum(rand(rng, F) * x^i for i in 0:t-1)
        is_irreducible(g) && return g
    end
end

# Binary vector of the m*t coordinates of the coefficients of a polynomial of
# degree less than t, as used for syndromes.
function _goppa_bits(f::FqPolyRingElem, t::Int)
    return reduce(vcat, [_gf2m_bits(coeff(f, j)) for j in 0:t-1])
end

# Polynomial of degree less than t with the coefficients given by `bits`.
function _goppa_poly(x::FqPolyRingElem, bits::AbstractVector{Int}, m::Int, t::Int)
    F = base_ring(parent(x))
    a = gen(F)
    return sum(_gf2m_element(F, a, bits[j*m+1:(j+1)*m]) * x^j for j in 0:t-1)
end

# Random binary Goppa code of length `n` over GF(2^m) correcting `t` errors.
# Codes whose parity check matrix does not have full rank m*t are rejected,
# so that the dimension is always n - m*t.
function _goppa_code(rng::AbstractRNG, m::Int, n::Int, t::Int)
    F, a = finite_field(2, m, :a)
    _, x = polynomial_ring(F, :x; cached = false)
    while true
        g = _random_irreducible(rng, x, t)
        idx = shuffle(rng, 0:2^m-1)[1:n]
        support = [_gf2m_element(F, a, [(i >> j) & 1 for j in 0:m-1]) for i in idx]
        cols = [_goppa_bits(invmod(x - s, g), t) for s in support]
        H = matrix(GF(2), reduce(hcat, cols))
        rank(H) == m * t || continue
        G = transpose(kernel(H; side = :right))
        return BinaryGoppaCode(m, t, support, g, H, G)
    end
end

# Square root of `f` modulo the irreducible polynomial `g` of degree t over
# GF(2^m): squaring is a bijection of a field with 2^(mt) elements, whose
# inverse is raising to the power 2^(mt-1).
_sqrt_mod(f::FqPolyRingElem, g::FqPolyRingElem, m::Int, t::Int) =
    powermod(f, ZZ(2)^(m * t - 1), g)

# Positions (between 1 and n) of the errors of a word with the binary syndrome
# `s` with respect to the parity check matrix of the code `C`. Patterson's
# algorithm computes the error locator sigma = prod (x - a_i) over the error
# positions i from the syndrome polynomial S = sigma' / sigma mod g. Returns
# `nothing` if there is no error vector of weight at most t with syndrome `s`.
function _goppa_decode(C::BinaryGoppaCode, s::AbstractVector{Int})
    m, t, g = C.m, C.t, C.g
    x = gen(parent(g))
    S = _goppa_poly(x, s, m, t)
    is_zero(S) && return Int[]
    T = invmod(S, g)
    if T == x
        sigma = x
    else
        # Find a, b with a = b * sqrt(T + x) mod g and deg(a) <= t/2,
        # deg(b) <= (t-1)/2 with the extended Euclidean algorithm.
        r0, r1 = g, _sqrt_mod(T + x, g, m, t)
        b0, b1 = zero(g), one(g)
        while degree(r1) > t ÷ 2
            q, r = divrem(r0, r1)
            r0, r1 = r1, r
            b0, b1 = b1, b0 - q * b1
        end
        sigma = r1^2 + x * b1^2
    end
    pos = [i for (i, a) in enumerate(C.support) if is_zero(evaluate(sigma, a))]
    return length(pos) == degree(sigma) ? pos : nothing
end
