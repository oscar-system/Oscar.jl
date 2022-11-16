@doc Markdown.doc"""
    domain(tm::ToricMorphism)

Return the domain of the toric morphism `tm`.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> domain(ToricIdentityMorphism(F4))
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor
```
"""
domain(tm::ToricMorphism) = tm.domain
export domain


@doc Markdown.doc"""
    codomain(tm::ToricMorphism)

Return the codomain of the toric morphism `tm`.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> codomain(ToricIdentityMorphism(F4))
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor
```
"""
codomain(tm::ToricMorphism) = tm.codomain
export codomain


@doc Markdown.doc"""
    image(tm::ToricMorphism)

Return the image of the toric morphism `tm`.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> image(ToricIdentityMorphism(F4))
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor
```
"""
image(tm::ToricMorphism) = tm.image
export image


@doc Markdown.doc"""
    grid_morphism(tm::ToricMorphism)

Return the underlying grid morphism of the toric morphism `tm`.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> grid_morphism(ToricIdentityMorphism(F4))
Map with following data
Domain:
=======
Abelian group with structure: Z^2
Codomain:
=========
Abelian group with structure: Z^2
```
"""
grid_morphism(tm::ToricMorphism) = tm.grid_morphism
export grid_morphism


@doc Markdown.doc"""
    morphism_on_torusinvariant_weil_divisor_group(tm::ToricMorphism)

For a given toric morphism `tm`, this method computes the corresponding
map of the torusinvariant Weil divisors.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> morphism_on_torusinvariant_weil_divisor_group(ToricIdentityMorphism(F4))
Map with following data
Domain:
=======
Abelian group with structure: Z^4
Codomain:
=========
Abelian group with structure: Z^4
```
"""
@attr GrpAbFinGenMap function morphism_on_torusinvariant_weil_divisor_group(tm::ToricMorphism)
    d = domain(tm)
    cod = codomain(tm)
    cod_rays = matrix(ZZ, rays(cod))
    images = matrix(ZZ, rays(d)) * matrix(grid_morphism(tm))
    mapping_matrix = matrix(ZZ, zeros(ZZ, rank(torusinvariant_weil_divisor_group(d)), 0))
    for i in 1:nrows(images)
      j = findfirst(x -> x == true, [contains(maximal_cones(cod)[j], [images[i,k] for k in 1:ncols(images)]) for j in 1:n_maximal_cones(cod)])
      m = vcat([Int(ray_indices(maximal_cones(cod))[j, k]) * cod_rays[k, :] for k in 1:nrays(cod)])
      mapping_matrix = hcat(mapping_matrix, solve(transpose(m), transpose(images[i, :])))
    end
    return hom(torusinvariant_weil_divisor_group(d), torusinvariant_weil_divisor_group(cod), mapping_matrix)
end
export morphism_on_torusinvariant_weil_divisor_group


@doc Markdown.doc"""
    morphism_on_cartier_divisor_group(tm::ToricMorphism)

For a given toric morphism `tm`, this method computes the corresponding
map of the Cartier divisors.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> morphism_on_cartier_divisor_group(ToricIdentityMorphism(F4))
Map with following data
Domain:
=======
Abelian group with structure: Z^4
Codomain:
=========
Abelian group with structure: Z^4
```
"""
@attr GrpAbFinGenMap function morphism_on_torusinvariant_cartier_divisor_group(tm::ToricMorphism)
    domain_variety = domain(tm)
    codomain_variety = codomain(tm)
    source_embedding = map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(domain_variety)
    morphism_of_weil_divisors = morphism_on_torusinvariant_weil_divisor_group(tm)
    return restrict_codomain(source_embedding * morphism_of_weil_divisors)
end
export morphism_on_torusinvariant_cartier_divisor_group
