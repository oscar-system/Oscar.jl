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

struct PuiseuxMPolyRing{T <: FieldElement} <: Ring
    underlyingPolynomialRing::MPolyRing

    function PuiseuxMPolyRing(R::MPolyRing)
        return new{elem_type(base_ring_type(R))}(R)
    end
end

mutable struct PuiseuxMPolyRingElem{T <: FieldElement} <: RingElem
    parent::PuiseuxMPolyRing{T}
    poly::MPolyRingElem
    shift::Vector{ZZRingElem}
    scale::ZZRingElem

    function PuiseuxMPolyRingElem(
        Kt::PuiseuxMPolyRing,
        f::MPolyRingElem,
        k::Vector{ZZRingElem} = zeros(ZZRingElem, nvars(Kt)),
        d::ZZRingElem = one(ZZ)
        )

        @req parent(f) == underlying_polynomial_ring(Kt) "polynomial must be in the underlying polynomial ring"
        @req d > 0 "scale must be positive"
        return new{elem_type(base_ring_type(parent(f)))}(Kt, f, k, d)
    end
end

function puiseux_polynomial_ring_elem(
    Kt::PuiseuxMPolyRing,
    f::MPolyRingElem,
    k::Vector{ZZRingElem} = zeros(ZZRingElem, nvars(Kt)),
    d::ZZRingElem = one(ZZ);
    skip_normalization::Bool = false
    )
    pf = PuiseuxMPolyRingElem(Kt, f, k, d)
    if !skip_normalization
        normalize!(pf)
    end
    return pf
end

#################################################################################
#
# Constructors
#
#################################################################################

function puiseux_polynomial_ring(K::Field, variableName::Vector{String})
    @req !isempty(variableName) "list of variables must not be empty"
    base_ring, _ = polynomial_ring(K, variableName)
    Kt = PuiseuxMPolyRing(base_ring)
    return Kt, gens(Kt)
end


#################################################################################
#
# Getters
#
#################################################################################

underlying_polynomial_ring(R::PuiseuxMPolyRing{T}) where T = R.underlyingPolynomialRing::mpoly_ring_type(T)
base_ring(R::PuiseuxMPolyRing) = base_ring(underlying_polynomial_ring(R))
coefficient_ring(R::PuiseuxMPolyRing) = base_ring(R)

Base.parent(f::PuiseuxMPolyRingElem) = f.parent
poly(f::PuiseuxMPolyRingElem{T}) where {T} = f.poly::mpoly_type(T)
shift(f::PuiseuxMPolyRingElem) = f.shift
scale(f::PuiseuxMPolyRingElem) = f.scale


#################################################################################
#
# Setters
#
#################################################################################

# WARNING: input is not assumed to be normalized
function normalize!(f::PuiseuxMPolyRingElem)
    if iszero(f)
        f.shift .= zeros(ZZRingElem, nvars(parent(f)))
        f.scale = one(ZZ)
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
function rescale(f::PuiseuxMPolyRingElem, newScale::ZZRingElem)
    @req newScale > 0 "new scale must be positive"

    # we assume f is normalized
    if newScale == scale(f)
        return f
    end

    newScaleMultipleOfCurrentScale, scaleQuotient = divides(newScale, scale(f))
    @req newScaleMultipleOfCurrentScale "new scale must be a multiple of the current scale"

    t = gens(underlying_polynomial_ring(parent(f)))
    newPoly = evaluate(poly(f), t .^ scaleQuotient)
    newShift = shift(f) * scaleQuotient
    return PuiseuxMPolyRingElem(parent(f), newPoly, newShift, newScale)
end


#################################################################################
#
# Conversions
#
#################################################################################

# The next function is required but not tested in AbstractAlgebra.
# The following code errors without it:
#   QQt,(t,) = puiseux_polynomial_ring(QQ,["t"]);
#   QQtx, x = polynomial_ring(QQt,3);
#   prod(x .^ rand(1:9,3)) # calls mul_johnson in AA/src/generic/MPoly,jl
function (R::PuiseuxMPolyRing)()
    return zero(R)
end

function (Kt::PuiseuxMPolyRing)(c::RingElement)
    return PuiseuxMPolyRingElem(Kt,underlying_polynomial_ring(Kt)(c))
end

function (Kt::PuiseuxMPolyRing{T})(ct::PuiseuxMPolyRingElem{T}) where T <: FieldElement
    return ct
end



#################################################################################
#
# Properties
#
#################################################################################

elem_type(::Type{PuiseuxMPolyRing{T}}) where T <: FieldElement = PuiseuxMPolyRingElem{T}
parent_type(::Type{PuiseuxMPolyRingElem{T}}) where T <: FieldElement = PuiseuxMPolyRing{T}
base_ring_type(::Type{PuiseuxMPolyRing{T}}) where T <: FieldElement = parent_type(T)

# The next function is required but not tested in AbstractAlgebra.
# The following code errors without it:
# K = algebraic_closure(QQ);
# Kz, z = polynomial_ring(K, "z");
# C = roots(rand(Int8)*z^2+rand(Int8)*z+rand(Int8))
# for _ in 1:99
#     C = vcat(C,roots(rand(Int8)*z^2+rand(Int8)*z+rand(Int8)))
# end
# Kt,(t,) = puiseux_polynomial_ring(K,["t"]);
# Ct = [ Kt(c) * t^rand(Int8) for c in C ]
# Ktx,x = polynomial_ring(Kt,3);
# f = sum([ rand(Ct) * prod(x .^ rand(1:9,3))  for _ in 1:9])
# evaluate(f, [Ktx(1), Ktx(1), Ktx(1)+x[3]]) # runs ^(::MPoly,Int) in AA/src/generic/MPoly,jl
characteristic(R::PuiseuxMPolyRing) = characteristic(coefficient_ring(R))

function Base.hash(f::PuiseuxMPolyRingElem, h::UInt)
    normalize!(f)
    return hash((parent(f), poly(f), shift(f), scale(f)), h)
end

gens(R::PuiseuxMPolyRing) = puiseux_polynomial_ring_elem.(Ref(R), gens(underlying_polynomial_ring(R)))
ngens(R::PuiseuxMPolyRing) = ngens(underlying_polynomial_ring(R))
nvars(R::PuiseuxMPolyRing) = nvars(underlying_polynomial_ring(R))
zero(R::PuiseuxMPolyRing) = puiseux_polynomial_ring_elem(R, zero(underlying_polynomial_ring(R)); skip_normalization=true)
one(R::PuiseuxMPolyRing) = puiseux_polynomial_ring_elem(R, one(underlying_polynomial_ring(R)); skip_normalization=true)
iszero(f::PuiseuxMPolyRingElem) = iszero(poly(f))
isone(f::PuiseuxMPolyRingElem) = isone(poly(f)) && iszero(shift(f)) && scale(f) == 1

function Base.:(==)(f::PuiseuxMPolyRingElem, g::PuiseuxMPolyRingElem)
    check_parent(f, g)
    return poly(f) == poly(g) && shift(f) == shift(g) && scale(f) == scale(g)
end

function Base.deepcopy_internal(f::PuiseuxMPolyRingElem, dict::IdDict)
    return puiseux_polynomial_ring_elem(parent(f), deepcopy_internal(poly(f), dict), deepcopy_internal(shift(f), dict), deepcopy_internal(scale(f), dict); skip_normalization=true)
end

coefficients(f::PuiseuxMPolyRingElem) = coefficients(poly(f))
exponents(f::PuiseuxMPolyRingElem) = [ (e + shift(f)) // scale(f) for e in exponents(poly(f)) ]
monomials(f::PuiseuxMPolyRingElem) = puiseux_polynomial_ring_elem.(Ref(parent(f)), monomials(poly(f)), Ref(shift(f)), Ref(scale(f)))

Base.length(f::PuiseuxMPolyRingElem) = length(poly(f))

function valuation(f::PuiseuxMPolyRingElem)
    @req nvars(parent(f)) == 1 "valuation is only defined for univariate Puiseux polynomials"
    if iszero(f)
        return PosInf()
    end
    return first(shift(f)) // scale(f)
end

is_univariate(R::PuiseuxMPolyRing) = is_univariate(underlying_polynomial_ring(R))
is_gen(f::PuiseuxMPolyRingElem) = is_gen(poly(f)) && iszero(shift(f)) && scale(f) == 1
is_monomial(f::PuiseuxMPolyRingElem) = is_monomial(poly(f))
is_unit(f::PuiseuxMPolyRingElem) = is_monomial(f) && is_unit(leading_coefficient(poly(f)))


#################################################################################
#
# Printing
#
#################################################################################

function Base.show(io::IO, R::PuiseuxMPolyRing)
    io = pretty(io)
    print(io, LowercaseOff(), "Puiseux polynomial ring over ", Lowercase(), coefficient_ring(R))
end

function Base.show(io::IO, f::PuiseuxMPolyRingElem)
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

function Base.:+(f::PuiseuxMPolyRingElem, g::PuiseuxMPolyRingElem)
    check_parent(f, g)
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
    fplusg = PuiseuxMPolyRingElem(
        parent(f),
        newPoly,
        newShift,
        newScale)
    normalize!(fplusg)
    return fplusg
end

function Base.:-(f::PuiseuxMPolyRingElem)
    return PuiseuxMPolyRingElem(parent(f), -poly(f), shift(f), scale(f))
end

function Base.:-(f::PuiseuxMPolyRingElem, g::PuiseuxMPolyRingElem)
    return f + (-g)
end

function Base.:*(f::PuiseuxMPolyRingElem, g::PuiseuxMPolyRingElem)
    check_parent(f, g)
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

function Base.:^(f::PuiseuxMPolyRingElem, a::QQFieldElem)

    if denominator(a) == 1
        return f^numerator(a)
    end

    @req length(f) == 1 "only monomials can be exponentiated to rational powers"
    @req isone(first(coefficients(f))) "only monomials with coefficient 1 can be exponentiated to rational powers"

    return puiseux_polynomial_ring_elem(
        parent(f),
        poly(f),
        shift(f)*numerator(a),
        scale(f)*denominator(a)
    )
end

function Base.:^(f::PuiseuxMPolyRingElem, a::Rational{Int})
    return f^(QQ(a))
end

function Base.:^(f::PuiseuxMPolyRingElem, a::ZZRingElem)
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

function Base.:^(f::PuiseuxMPolyRingElem, a::Integer)
    return f^(ZZ(a))
end

# TODO: these // should be changed to divexact


function Base.:(//)(f::PuiseuxMPolyRingElem{K}, a::K) where K <: FieldElement
    @req !iszero(a) "division by zero"
    return puiseux_polynomial_ring_elem(parent(f), poly(f)*1//a, shift(f), scale(f); skip_normalization=true)
end

function Base.:(//)(f::PuiseuxMPolyRingElem, a::Int)
    # TODO: delegate to divexact
    return f//coefficient_ring(f)(a)
end

# The next function is required but not tested in AbstractAlgebra.
# The following code errors without it (after implenenting characteristic above):
# K = algebraic_closure(QQ);
# Kz, z = polynomial_ring(K, "z");
# C = roots(rand(Int8)*z^2+rand(Int8)*z+rand(Int8))
# for _ in 1:99
#     C = vcat(C,roots(rand(Int8)*z^2+rand(Int8)*z+rand(Int8)))
# end
# Kt,(t,) = puiseux_polynomial_ring(K,["t"]);
# Ct = [ Kt(c) * t^rand(Int8) for c in C ]
# Ktx,x = polynomial_ring(Kt,3);
# f = sum([ rand(Ct) * prod(x .^ rand(1:9,3))  for _ in 1:9])
# evaluate(f, [Ktx(1), Ktx(1), Ktx(1)+x[3]]) # runs pow_fps in AA/src/generic/MPoly.jl
function divexact(f::PuiseuxMPolyRingElem, g::PuiseuxMPolyRingElem)
    check_parent(f, g)
    @req !iszero(g) "division by zero"

    if iszero(f)
        return zero(parent(f))
    end
    if isone(g)
        return f
    end

    @assert length(g) == 1 "can only divide by monomials"

    # subtract shifts, multiply scales and divide poly(f) by coefficient of poly(g)
    newShift = shift(f)*scale(g) - shift(g)*scale(f)
    newScale = scale(f)*scale(g)
    xPrimes = gens(parent(poly(f)))
    # TODO: the '/' in the line below should be divexact of polys
    newPoly = evaluate(poly(f), xPrimes .^ scale(g)) / first(coefficients(poly(g)))

    return puiseux_polynomial_ring_elem(parent(f), newPoly, newShift, newScale)
end

# The following function is required for running the Conformance Tests
function ConformanceTests.generate_element(R::PuiseuxMPolyRing)
    f = ConformanceTests.generate_element(underlying_polynomial_ring(R))
    shift = [rand(ZZ, -10:10) for i in 1:nvars(R)]
    scale = rand(ZZ, 1:10)
    return puiseux_polynomial_ring_elem(R, f, shift, scale)
end
