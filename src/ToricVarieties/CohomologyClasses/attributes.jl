########################
# Attributes of cohomology classes
########################

@doc Markdown.doc"""
    toric_variety(c::CohomologyClass)

Return the normal toric variety of the cohomology class `c`.
"""
function toric_variety(c::CohomologyClass)
    return c.toric_variety
end
export toric_variety


@doc Markdown.doc"""
    coefficients(c::CohomologyClass)

Return the coefficients of the cohomology class `c`.
"""
function coefficients(c::CohomologyClass)
    return [coefficient_ring(toric_variety(c))(k) for k in c.coeffs]
end
export coefficients


@doc Markdown.doc"""
    exponents(c::CohomologyClass)

Return the exponents of the cohomology class `c`.
"""
function exponents(c::CohomologyClass)
    return c.exponents
end
export exponents


@doc Markdown.doc"""
    polynomial(c::CohomologyClass)

Return the polynomial in the cohomology ring of the normal
toric variety `toric_variety(c)` which corresponds to `c`.
"""
function polynomial(c::CohomologyClass)
    return polynomial(c, cohomology_ring(toric_variety(c)))
end
export polynomial


@doc Markdown.doc"""
    polynomial(c::CohomologyClass, ring::MPolyQuo)

Returns the polynomial in `ring` corresponding 
to the cohomology class `c`.
"""
function polynomial(c::CohomologyClass, ring::MPolyQuo)
    coeffs = coefficients(c)
    if length(coeffs) == 0
        return zero(ring)
    end
    expos = exponents(c)
    indets = gens(ring)
    monoms = [prod(indets[j]^expos[k,j] for j in 1:ncols(expos)) for k in 1:nrows(expos)]
    return sum(coeffs[k]*monoms[k] for k in 1:length(monoms))
end
export polynomial
