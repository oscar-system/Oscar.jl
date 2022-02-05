#####################
# 1. Defining data of a line bundle
#####################

@doc Markdown.doc"""
    divisor_class(l::ToricLineBundle)

Return the divisor class which defines the toric line bundle `l`.
"""
function divisor_class(l::ToricLineBundle)
    return l.divisor_class
end
export divisor_class


@doc Markdown.doc"""
    variety(l::ToricLineBundle)

Return the toric variety over which the toric line bundle `l` is defined.
"""
function variety(l::ToricLineBundle)
    return l.variety
end
export variety


@doc Markdown.doc"""
    toric_divisor(l::ToricLineBundle)

Returns a divisor corresponding to the toric line bundle `l`.
"""
@attr ToricDivisor function toric_divisor(l::ToricLineBundle)
    class = divisor_class(l)
    map1 = map_from_cartier_divisor_group_to_picard_group(variety(l))
    map2 = map = map_from_cartier_divisor_group_to_torus_invariant_divisor_group(variety(l))
    image = map2(preimage(map1, class)).coeff
    coeffs = vec([fmpz(x) for x in image])
    td = ToricDivisor(variety(l), coeffs)
    set_attribute!(td, :iscartier, true)
    return td
end
export toric_divisor


@doc Markdown.doc"""
    degree(l::ToricLineBundle)

Returns the degree of the toric line bundle `l`.
"""
@attr fmpz function degree(l::ToricLineBundle)
    return sum(coefficients(toric_divisor(l)))
end
export degree
