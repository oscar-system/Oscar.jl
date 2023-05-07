```@meta
CurrentModule = Oscar
```

# Normal Toric Schemes

## Constructors

We can construct a toric scheme as follows:
```julia
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> toric_scheme = ToricCoveredScheme(P2)
Scheme of a toric variety with fan spanned by RayVector{QQFieldElem}[[1, 0], [0, 1], [-1, -1]]
```

## Special constructors

We support constructors for a couple of famous geometries:
```@docs
affine_space(::Type{ToricCoveredScheme}, d::Int; set_attributes::Bool = true)
projective_space(::Type{ToricCoveredScheme}, d::Int; set_attributes::Bool = true)
weighted_projective_space(::Type{ToricCoveredScheme}, w::Vector{T}; set_attributes::Bool = true) where {T <: IntegerUnion}
hirzebruch_surface(::Type{ToricCoveredScheme}, r::Int; set_attributes::Bool = true)
del_pezzo_surface(::Type{ToricCoveredScheme}, b::Int; set_attributes::Bool = true)
```


## Attributes

We overload the following attributes from the underlying toric variety:
* ``affine_open_covering(X::ToricCoveredScheme)``,
* ``character_lattice(X::ToricCoveredScheme)``,
* ``character_to_rational_function(R::MPolyRing, X::ToricCoveredScheme, character::Vector{ZZRingElem})``,
* ``character_to_rational_function(X::ToricCoveredScheme, character::Vector{ZZRingElem})``,
* ``class_group(X::ToricCoveredScheme)``,
* ``coefficient_ring(X::ToricCoveredScheme)``,
* ``coordinate_names(X::ToricCoveredScheme)``,
* ``coordinate_names_of_torus(X::ToricCoveredScheme)``,
* ``coordinate_ring_of_torus(R::MPolyRing, X::ToricCoveredScheme)``,
* ``coordinate_ring_of_torus(X::ToricCoveredScheme)``,
* ``cox_ring(R::MPolyRing, X::ToricCoveredScheme)``,
* ``cox_ring(X::ToricCoveredScheme)``,
* ``dim(X::ToricCoveredScheme)``,
* ``dim_of_torusfactor(X::ToricCoveredScheme)``,
* ``euler_characteristic(X::ToricCoveredScheme)``,
* ``fan(X::ToricCoveredScheme)``,
* ``ideal_of_linear_relations(R::MPolyRing, X::ToricCoveredScheme)``,
* ``ideal_of_linear_relations(X::ToricCoveredScheme)``,
* ``irrelevant_ideal(R::MPolyRing, X::ToricCoveredScheme)``,
* ``irrelevant_ideal(X::ToricCoveredScheme)``,
* ``map_from_character_lattice_to_torusinvariant_weil_divisor_group(X::ToricCoveredScheme)``,
* ``map_from_torusinvariant_cartier_divisor_group_to_picard_group(X::ToricCoveredScheme)``,
* ``map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(X::ToricCoveredScheme)``,
* ``map_from_torusinvariant_weil_divisor_group_to_class_group(X::ToricCoveredScheme)``,
* ``mori_cone(X::ToricCoveredScheme)``,
* ``nef_cone(X::ToricCoveredScheme)``,
* ``picard_group(X::ToricCoveredScheme)``,
* ``set_coordinate_names(X::ToricCoveredScheme)``,
* ``set_coordinate_names_of_torus(X::ToricCoveredScheme)``,
* ``stanley_reisner_ideal(R::MPolyRing, X::ToricCoveredScheme)``,
* ``stanley_reisner_ideal(X::ToricCoveredScheme)``,
* ``torusinvariant_cartier_divisor_group(X::ToricCoveredScheme)``,
* ``torusinvariant_prime_divisors(X::ToricCoveredScheme)``,
* ``torusinvariant_weil_divisor_group(X::ToricCoveredScheme)``,
* ``underlying_toric_variety(X::ToricCoveredScheme)``.
It is possible to forget the toric origin, i.e. one can
turn a toric scheme into a scheme. This is achieved
with the following method:
```@docs
underlying_scheme(X::ToricCoveredScheme)
```


## Properties

The following properties are overloaded from the underlying toric variety:
* ``is_finalized(X::ToricCoveredScheme)``,
* ``is_normal(X::ToricCoveredScheme)``,
* ``is_affine(X::ToricCoveredScheme)``,
* ``is_projective(X::ToricCoveredScheme)``,
* ``is_projective_space(X::ToricCoveredScheme)``,
* ``is_smooth(X::ToricCoveredScheme)``,
* ``is_complete(X::ToricCoveredScheme)``,
* ``has_torusfactor(X::ToricCoveredScheme)``,
* ``is_orbifold(X::ToricCoveredScheme)``,
* ``is_simplicial(X::ToricCoveredScheme)``,
* ``is_gorenstein(X::ToricCoveredScheme)``,
* ``is_q_gorenstein(X::ToricCoveredScheme)``,
* ``is_fano(X::ToricCoveredScheme)``.
