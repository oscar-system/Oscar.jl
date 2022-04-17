#####################
# 1. Defining data of a line bundle
#####################

@doc Markdown.doc"""
    divisor_class(l::ToricLineBundle)

Return the divisor class which defines the toric line bundle `l`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> l = ToricLineBundle(v, [fmpz(2)])
A toric line bundle on a normal toric variety

julia> divisor_class(l)
Element of
GrpAb: Z
with components [2]
```
"""
function divisor_class(l::ToricLineBundle)
    return l.divisor_class
end
export divisor_class


@doc Markdown.doc"""
    toric_variety(l::ToricLineBundle)

Return the toric variety over which the toric line bundle `l` is defined.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> l = ToricLineBundle(v, [fmpz(2)])
A toric line bundle on a normal toric variety

julia> toric_variety(l)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor
```
"""
function toric_variety(l::ToricLineBundle)
    return l.toric_variety
end
export toric_variety


@doc Markdown.doc"""
    toric_divisor(l::ToricLineBundle)

Return a divisor corresponding to the toric line bundle `l`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> l = ToricLineBundle(v, [fmpz(2)])
A toric line bundle on a normal toric variety

julia> toric_divisor(l)
A torus-invariant, cartier, non-prime divisor on a normal toric variety

julia> iscartier(toric_divisor(l))
true
```
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

Return the degree of the toric line bundle `l`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> l = ToricLineBundle(v, [fmpz(2)])
A toric line bundle on a normal toric variety

julia> degree(l)
2
```
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

Return a basis of the global sections of the toric line bundle `l` in terms of rational functions.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> l = ToricLineBundle(v, [fmpz(2)])
A toric line bundle on a normal toric variety

julia> basis_of_global_sections_via_rational_functions(l)
6-element Vector{MPolyQuoElem{fmpq_mpoly}}:
 x1_^2
 x2*x1_^2
 x2^2*x1_^2
 x1_
 x2*x1_
 1
```
"""
@attr function basis_of_global_sections_via_rational_functions(l::ToricLineBundle)
    characters = matrix(ZZ, lattice_points(polyhedron(toric_divisor(l))))
    return [character_to_rational_function(toric_variety(l), vec([fmpz(c) for c in characters[i,:]])) for i in 1:nrows(characters)]
end
export basis_of_global_sections_via_rational_functions


@doc Markdown.doc"""
    basis_of_global_sections_via_homogeneous_component(l::ToricLineBundle)

Return a basis of the global sections of the toric line bundle `l`
in terms of a homogeneous component of the Cox ring of `toric_variety(l)`.
For convenience, this method can also be called via
`basis_of_global_sections(l::ToricLineBundle)`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> l = ToricLineBundle(v, [fmpz(2)])
A toric line bundle on a normal toric variety

julia> basis_of_global_sections_via_homogeneous_component(l)
6-element Vector{MPolyElem_dec{fmpq, fmpq_mpoly}}:
 x3^2
 x2*x3
 x2^2
 x1*x3
 x1*x2
 x1^2

julia> basis_of_global_sections(l)
6-element Vector{MPolyElem_dec{fmpq, fmpq_mpoly}}:
 x3^2
 x2*x3
 x2^2
 x1*x3
 x1*x2
 x1^2
```
"""
@attr function basis_of_global_sections_via_homogeneous_component(l::ToricLineBundle)
    hc = homogeneous_component(cox_ring(toric_variety(l)), divisor_class(l))
    generators = gens(hc[1])
    return [hc[2](x) for x in generators]
end
basis_of_global_sections(l::ToricLineBundle) = basis_of_global_sections_via_homogeneous_component(l)
export basis_of_global_sections_via_homogeneous_component, basis_of_global_sections
