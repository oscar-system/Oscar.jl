#################################################################################
#
# Puiseux polynomials
# ===================
#
# A multivariate Puiseux polynomial ring is just a wrapper around a normal
# multivariate polynomial ring for the normal poly below.
#
# A Puiseux polynomial can be written of the form
#
#   t^{k//d} * (c_0 * t^(0//d) + ... + c_r * t^(r//d)).
#
# We store it as
# - a normal poly: c_0 * t^0 + ... + c_r * t^r
# - a shift: k
# - a scale: d
#
# The representation above is normalized if
# - gcd(d, exponents of poly) = 1
# - c_0 != 0
# Every Puiseux polynomial has a unique normalized representation.
#
# With the exception of normalized! and rescale, all function inputs are assumed
# to be normalized and all function outputs will be normalized.
#
#################################################################################

struct MPuiseuxPolyRing{T <: FieldElement} <: Ring
    underlyingPolynomialRing::MPolyRing

    function MPuiseuxPolyRing(R::MPolyRing)
        return new{elem_type(base_ring_type(R))}(R)
    end
end

mutable struct MPuiseuxPolyRingElem{T <: FieldElement} <: RingElem
    parent::MPuiseuxPolyRing{T}
    poly::MPolyRingElem
    shift::Vector{ZZRingElem}
    scale::ZZRingElem

    function MPuiseuxPolyRingElem(
        Kt::MPuiseuxPolyRing,
        f::MPolyRingElem,
        k::Vector{ZZRingElem} = zeros(ZZRingElem, nvars(Kt)),
        d::ZZRingElem = one(ZZ)
        )

        @assert parent(f) == underlying_polynomial_ring(Kt) "polynomial must be in the underlying polynomial ring"
        @assert d > 0 "scale must be positive"
        return new{elem_type(base_ring_type(parent(f)))}(Kt, f, k, d)
    end
end

function puiseux_polynomial_ring_elem(
    Kt::MPuiseuxPolyRing,
    f::MPolyRingElem,
    k::Vector{ZZRingElem} = zeros(ZZRingElem, nvars(Kt)),
    d::ZZRingElem = one(ZZ);
    skip_normalization::Bool = false
    )
    pf = MPuiseuxPolyRingElem(Kt, f, k, d)
    if !skip_normalization
        normalize!(pf)
    end
    return pf
end

#################################################################################
#
# Constructors
#
##################################################################################

function puiseux_polynomial_ring(K::Field, variableName::Vector{String})
    @assert !isempty(variableName) "list of variables must not be empty"
    base_ring, _ = polynomial_ring(K, variableName)
    Kt = MPuiseuxPolyRing(base_ring)
    return Kt, gens(Kt)
end


#################################################################################
#
# Getters
#
#################################################################################

underlying_polynomial_ring(R::MPuiseuxPolyRing) = R.underlyingPolynomialRing
base_ring(R::MPuiseuxPolyRing) = base_ring(underlying_polynomial_ring(R))
coefficient_ring(R::MPuiseuxPolyRing) = base_ring(R)

Base.parent(f::MPuiseuxPolyRingElem) = f.parent
poly(f::MPuiseuxPolyRingElem) = f.poly
shift(f::MPuiseuxPolyRingElem) = f.shift
scale(f::MPuiseuxPolyRingElem) = f.scale


#################################################################################
#
# Setters
#
#################################################################################

# WARNING: input is not assumed to be normalized
function normalize!(f::MPuiseuxPolyRingElem)
    if iszero(f)
        return false
    end

    underlyingPolynomialRing = underlying_polynomial_ring(parent(f))

    # make sure shift is correct, i.e., for every variable poly(f) has a
    # monomial that is constant in that variable
    shiftDifference = [ reduce(min,[e[i] for e in exponents(poly(f))]) for i in 1:nvars(parent(poly(f))) ]
    if !iszero(shiftDifference)
        t = gens(underlyingPolynomialRing)
        f.poly = div(f.poly, prod(t.^(shiftDifference)))
        f.shift += shiftDifference
    end

    # make sure scale is correct,
    # i.e., gcd of numerators (= exponents of poly + shift) and denominatos (= scale) is 1
    gcdExponents = gcd(vcat(scale(f), reduce(vcat,[e+shift(f) for e in exponents(poly(f))])))
    if gcdExponents > 1
        f.poly = sum([c*monomial(underlyingPolynomialRing, Int.(div.(e, gcdExponents))) for (c, e) in zip(coefficients(poly(f)), exponents(poly(f)))])
        f.shift = div.(f.shift, gcdExponents)
        f.scale = div(f.scale, gcdExponents)
    end

    return !iszero(shiftDifference) || gcdExponents > 1
end

# WARNING: output may not be normalized
function rescale(f::MPuiseuxPolyRingElem, newScale::ZZRingElem)
    @assert newScale > 0 "new scale must be positive"

    # we assume f is normalized
    if newScale == scale(f)
        return f
    end

    newScaleMultipleOfCurrentScale, scaleQuotient = divides(newScale, scale(f))
    @assert newScaleMultipleOfCurrentScale "new scale must be a multiple of the current scale"

    t = gens(underlying_polynomial_ring(parent(f)))
    newPoly = evaluate(poly(f), t .^ scaleQuotient)
    newShift = shift(f) * scaleQuotient
    return MPuiseuxPolyRingElem(parent(f), newPoly, newShift, newScale)
end


#################################################################################
#
# Conversions
#
#################################################################################

function (Kt::MPuiseuxPolyRing)(c::Int)
    return MPuiseuxPolyRingElem(Kt,underlying_polynomial_ring(Kt)(c))
end

function (Kt::MPuiseuxPolyRing)(c::Rational{Int})
    return MPuiseuxPolyRingElem(Kt,underlying_polynomial_ring(Kt)(c))
end

function (Kt::MPuiseuxPolyRing)(c::RingElem)
    return MPuiseuxPolyRingElem(Kt,underlying_polynomial_ring(Kt)(c))
end

function (Kt::MPuiseuxPolyRing{T})(c::T) where T <: FieldElement
    return MPuiseuxPolyRingElem(Kt,underlying_polynomial_ring(Kt)(c))
end
function (Kt::MPuiseuxPolyRing{T})(ct::MPuiseuxPolyRingElem{T}) where T <: FieldElement
    return ct
end

#################################################################################
#
# Properties
#
#################################################################################

elem_type(::Type{MPuiseuxPolyRing{T}}) where T <: FieldElement = MPuiseuxPolyRingElem{T}
parent_type(::Type{MPuiseuxPolyRingElem{T}}) where T <: FieldElement = MPuiseuxPolyRing{T}
base_ring_type(::Type{MPuiseuxPolyRing{T}}) where T <: FieldElement = parent_type(T)

function Base.hash(f::MPuiseuxPolyRingElem, h::UInt)
    normalize!(f)
    return hash((parent(f), poly(f), shift(f), scale(f)), h)
end

gens(R::MPuiseuxPolyRing) = puiseux_polynomial_ring_elem.(Ref(R), gens(underlying_polynomial_ring(R)))
ngens(R::MPuiseuxPolyRing) = ngens(underlying_polynomial_ring(R))
nvars(R::MPuiseuxPolyRing) = nvars(underlying_polynomial_ring(R))
zero(R::MPuiseuxPolyRing) = puiseux_polynomial_ring_elem(R, zero(underlying_polynomial_ring(R)); skip_normalization=true)
one(R::MPuiseuxPolyRing) = puiseux_polynomial_ring_elem(R, one(underlying_polynomial_ring(R)); skip_normalization=true)
iszero(f::MPuiseuxPolyRingElem) = iszero(poly(f))
isone(f::MPuiseuxPolyRingElem) = isone(poly(f)) && shift(f) == 0 && scale(f) == 1

function Base.:(==)(f::MPuiseuxPolyRingElem, g::MPuiseuxPolyRingElem)
    @assert parent(f) == parent(g) "elements must be in the same ring"
    return poly(f) == poly(g) && shift(f) == shift(g) && scale(f) == scale(g)
end

coefficients(f::MPuiseuxPolyRingElem) = coefficients(poly(f))
exponents(f::MPuiseuxPolyRingElem) = [ (e + shift(f)) // scale(f) for e in exponents(poly(f)) ]
monomials(f::MPuiseuxPolyRingElem) = puiseux_polynomial_ring_elem.(Ref(parent(f)), monomials(poly(f)), Ref(shift(f)), Ref(scale(f)))

Base.length(f::MPuiseuxPolyRingElem) = length(poly(f))

function valuation(f::MPuiseuxPolyRingElem)
    @assert nvars(parent(f)) == 1 "valuation is only defined for univariate Puiseux polynomials"
    if iszero(f)
        return PosInf()
    end
    return first(shift(f)) // scale(f)
end


#################################################################################
#
# Printing
#
#################################################################################

function Base.show(io::IO, R::MPuiseuxPolyRing)
    print(io, "Puiseux polynomial ring over ", coefficient_ring(R))
end

function Base.show(io::IO, f::MPuiseuxPolyRingElem)
    if iszero(f)
        print(io, "0")
        return
    elseif isone(f)
        print(io, "1")
        return
    end

    t = gens(underlying_polynomial_ring(parent(f)))
    termStrings = []
    for (c, e) in zip(coefficients(poly(f)), exponents(poly(f)))
        exponentVector = (e + shift(f)) .// scale(f)
        if iszero(exponentVector)
            push!(termStrings, string(c))
        else
            monomialString = ""
            for (ti, ei) in zip(t, exponentVector)
                if iszero(ei)
                    continue
                elseif isone(ei)
                    monomialString *= "*" * string(ti)
                elseif denominator(ei) == 1 && numerator(ei) > 0
                    monomialString *= "*" * string(ti) * "^" * string(numerator(ei))
                else
                    monomialString *= "*" * string(ti) * "^(" * string(ei) * ")"
                end
            end
            if isone(c)
                push!(termStrings, monomialString[2:end]) # remove leading "*"
            else
                push!(termStrings, string(c) * monomialString)
            end
        end
    end
    print(io, join(termStrings, " + "))
end


#################################################################################
#
# Operations
#
#################################################################################

function Base.:+(f::MPuiseuxPolyRingElem, g::MPuiseuxPolyRingElem)
    @assert parent(f) == parent(g) "elements must be in the same ring"
    if iszero(f)
        return g
    elseif iszero(g)
        return f
    end

    # rescale both to the lcm of their scales
    newScale = lcm(scale(f), scale(g))
    frescaled = rescale(f, newScale)
    grescaled = rescale(g, newScale)

    # add the polynomials, adjusting for shifts
    newShift = min.(shift(frescaled), shift(grescaled))
    t = gens(underlying_polynomial_ring(parent(f)))
    newPoly = prod(t.^(shift(frescaled)-newShift))*poly(frescaled) + prod(t.^(shift(grescaled)-newShift))*poly(grescaled)

    # normalize output, in case of cancellations
    fplusg = MPuiseuxPolyRingElem(
        parent(f),
        newPoly,
        newShift,
        newScale)
    normalize!(fplusg)
    return fplusg
end

function Base.:-(f::MPuiseuxPolyRingElem)
    return MPuiseuxPolyRingElem(parent(f), -poly(f), shift(f), scale(f))
end

function Base.:-(f::MPuiseuxPolyRingElem, g::MPuiseuxPolyRingElem)
    return f + (-g)
end

function Base.:*(f::MPuiseuxPolyRingElem, g::MPuiseuxPolyRingElem)
    @assert parent(f) == parent(g) "elements must be in the same ring"
    if iszero(f) || iszero(g)
        return zero(parent(f))
    elseif isone(f)
        return g
    elseif isone(g)
        return f
    end

    # add shifts, multiply scales and polys
    newShift = shift(f)*scale(g) + shift(g)*scale(f)
    newScale = scale(f)*scale(g)
    xPrimes = gens(parent(poly(f)))
    newPoly = evaluate(poly(f),xPrimes .^ scale(g))*evaluate(poly(g),xPrimes .^ scale(f))

    return puiseux_polynomial_ring_elem(parent(f), newPoly, newShift, newScale)
end

function Base.:^(f::MPuiseuxPolyRingElem, a::QQFieldElem)

    if denominator(a) == 1
        return f^numerator(a)
    end

    @assert length(f) == 1 "only monomials can be exponentiated to rational powers"
    @assert isone(first(coefficients(f))) "only monomials with coefficient 1 can be exponentiated to rational powers"

    return puiseux_polynomial_ring_elem(
        parent(f),
        poly(f),
        shift(f)*numerator(a),
        scale(f)*denominator(a)
    )
end

function Base.:^(f::MPuiseuxPolyRingElem, a::Rational{Int})
    return f^(QQ(a))
end

function Base.:^(f::MPuiseuxPolyRingElem, a::ZZRingElem)
    if a == 0
        return one(parent(f))
    end
    if a == 1
        return f
    end

    if a < 0
        # test whether f is a monomial
        @assert length(f) == 1 "only monomials can be exponentiated to negative powers"
        R = parent(f)
        Runder = underlying_polynomial_ring(R)
        c = first(coefficients(f))
        return puiseux_polynomial_ring_elem(
            R,
            Runder(c^a),
            shift(f)*a,
            scale(f)
        )
    end

    return reduce(*, [ f for i in 1:a ])
end

function Base.:^(f::MPuiseuxPolyRingElem, a::Int)
    return f^(ZZ(a))
end

function Base.:(//)(f::MPuiseuxPolyRingElem{K}, a::K) where K <: FieldElement
    @assert !iszero(a) "division by zero"
    return puiseux_polynomial_ring_elem(parent(f), poly(f)*1//a, shift(f), scale(f); skip_normalization=true)
end

function Base.:(//)(f::MPuiseuxPolyRingElem, a::Int)
    return f//coefficient_ring(f)(a)
end
