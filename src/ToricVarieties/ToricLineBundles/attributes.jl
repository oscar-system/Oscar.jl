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
    toric_variety(l::ToricLineBundle)

Return the toric variety over which the toric line bundle `l` is defined.
"""
function toric_variety(l::ToricLineBundle)
    return l.toric_variety
end
export toric_variety


@doc Markdown.doc"""
    toric_divisor(l::ToricLineBundle)

Returns a divisor corresponding to the toric line bundle `l`.
"""
@attr ToricDivisor function toric_divisor(l::ToricLineBundle)
    class = divisor_class(l)
    map1 = map_from_torusinvariant_cartier_divisor_group_to_picard_group(toric_variety(l))
    map2 = map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(toric_variety(l))
    image = map2(preimage(map1, class)).coeff
    coeffs = vec([fmpz(x) for x in image])
    td = ToricDivisor(toric_variety(l), coeffs)
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


#####################
# 2. Basis of global sections
#####################

@doc Markdown.doc"""
    basis_of_global_sections_via_rational_functions(l::ToricLineBundle)

Returns a basis of the global sections of the toric line bundle `l`.
"""
@attr function basis_of_global_sections_via_rational_functions(l::ToricLineBundle)
    characters = matrix(ZZ, lattice_points(polyhedron(toric_divisor(l))))
    return [character_to_rational_function(toric_variety(l), vec([fmpz(c) for c in characters[i,:]])) for i in 1:nrows(characters)]
end
export basis_of_global_sections_via_rational_functions


@doc Markdown.doc"""
    basis_of_global_sections_via_homogeneous_component(l::ToricLineBundle)

Returns a basis of the global sections of the toric line bundle `l`,
which is computed from a homogeneous component of the cox ring.
"""
@attr function basis_of_global_sections_via_homogeneous_component(l::ToricLineBundle)
    hc = homogeneous_component(cox_ring(toric_variety(l)), divisor_class(l))
    generators = gens(hc[1])
    return [hc[2](x) for x in generators]
end
export basis_of_global_sections_via_homogeneous_component
