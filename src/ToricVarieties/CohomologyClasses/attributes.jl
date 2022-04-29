########################
# Attributes of cohomology classes
########################

@doc Markdown.doc"""
    toric_variety(c::CohomologyClass)

Return the normal toric variety of the cohomology class `c`.

# Examples
```jldoctest
julia> dP2 = del_pezzo(2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> d = ToricDivisor(dP2, [1,2,3,4,5])
A torus-invariant, non-prime divisor on a normal toric variety

julia> cc = CohomologyClass(d)
A cohomology class on a normal toric variety given by 6*x3 + e1 + 7*e2

julia> toric_variety(cc)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety over QQ without torusfactor
```
"""
toric_variety(c::CohomologyClass) = c.v
export toric_variety


@doc Markdown.doc"""
    coefficients(c::CohomologyClass)

Return the coefficients of the cohomology class `c`.

# Examples
```jldoctest
julia> dP2 = del_pezzo(2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> d = ToricDivisor(dP2, [1,2,3,4,5])
A torus-invariant, non-prime divisor on a normal toric variety

julia> cc = CohomologyClass(d)
A cohomology class on a normal toric variety given by 6*x3 + e1 + 7*e2

julia> coefficients(cc)
3-element Vector{fmpq}:
 6
 1
 7
```
"""
coefficients(c::CohomologyClass) = [coefficient_ring(toric_variety(c))(k) for k in coefficients(polynomial(c).f)]
export coefficients


@doc Markdown.doc"""
    exponents(c::CohomologyClass)

Return the exponents of the cohomology class `c`.

# Examples
```jldoctest
julia> dP2 = del_pezzo(2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> d = ToricDivisor(dP2, [1,2,3,4,5])
A torus-invariant, non-prime divisor on a normal toric variety

julia> cc = CohomologyClass(d)
A cohomology class on a normal toric variety given by 6*x3 + e1 + 7*e2

julia> exponents(cc)
[0   0   1   0   0]
[0   0   0   1   0]
[0   0   0   0   1]
```
"""
exponents(c::CohomologyClass) = matrix(ZZ,[k for k in exponent_vectors(polynomial(c).f)])
export exponents


@doc Markdown.doc"""
    polynomial(c::CohomologyClass)

Return the polynomial in the cohomology ring of the normal
toric variety `toric_variety(c)` which corresponds to `c`.

# Examples
```jldoctest
julia> dP2 = del_pezzo(2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> d = ToricDivisor(dP2, [1,2,3,4,5])
A torus-invariant, non-prime divisor on a normal toric variety

julia> cc = CohomologyClass(d)
A cohomology class on a normal toric variety given by 6*x3 + e1 + 7*e2

julia> polynomial(cc)
6*x3 + e1 + 7*e2
```
"""
polynomial(c::CohomologyClass) = c.p
export polynomial


@doc Markdown.doc"""
    polynomial(c::CohomologyClass, ring::MPolyQuo)

Returns the polynomial in `ring` corresponding
to the cohomology class `c`.

# Examples
```jldoctest
julia> dP2 = del_pezzo(2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> d = ToricDivisor(dP2, [1,2,3,4,5])
A torus-invariant, non-prime divisor on a normal toric variety

julia> cc = CohomologyClass(d)
A cohomology class on a normal toric variety given by 6*x3 + e1 + 7*e2

julia> R, _ = PolynomialRing(QQ, 5)
(Multivariate Polynomial Ring in x1, x2, x3, x4, x5 over Rational Field, fmpq_mpoly[x1, x2, x3, x4, x5])

julia> (x1,x2,x3,x4,x5) = gens(R)
5-element Vector{fmpq_mpoly}:
 x1
 x2
 x3
 x4
 x5

julia> sr_and_linear_relation_ideal = ideal([x1*x3, x1*x5, x2*x4, x2*x5, x3*x4, x1 + x2 - x5, x2 + x3 - x4 - x5])
ideal(x1*x3, x1*x5, x2*x4, x2*x5, x3*x4, x1 + x2 - x5, x2 + x3 - x4 - x5)

julia> R_quo = quo(R, sr_and_linear_relation_ideal)[1]
Quotient of Multivariate Polynomial Ring in x1, x2, x3, x4, x5 over Rational Field by ideal(x1*x3, x1*x5, x2*x4, x2*x5, x3*x4, x1 + x2 - x5, x2 + x3 - x4 - x5)

julia> polynomial(cc, R_quo)
6*x3 + x4 + 7*x5
```
"""
function polynomial(c::CohomologyClass, ring::MPolyQuo)
    p = polynomial(c)
    if iszero(p)
        return zero(ring)
    end
    coeffs = [k for k in coefficients(p.f)]
    expos = matrix(ZZ,[k for k in exponent_vectors(p.f)])
    indets = gens(ring)
    monoms = [prod(indets[j]^expos[k,j] for j in 1:ncols(expos)) for k in 1:nrows(expos)]
    return sum(coeffs[k]*monoms[k] for k in 1:length(monoms))
end
export polynomial
