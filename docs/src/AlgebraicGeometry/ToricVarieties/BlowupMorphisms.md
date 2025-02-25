```@meta
CurrentModule = Oscar
```

# Toric Blowups (Experimental)

It is a common goal in algebraic geometry to resolve singularities.
Certainly, toric varieties and their subvarieties are no exception and
we provide a growing set of functionality for such tasks.

We focus on toric blowups given by a star subdivision of a
polyhedral fan along a primitive vector, see 11.1 Star Subdivisions in
[CLS11](@cite).

Given a homogeneous ideal $I$ in the Cox ring of $X$, it is also possible to
reroute to covered schemes and construct the blowup by
`blow_up(ideal_sheaf(X, I))`.
Strict transforms (and similarly total transforms) of ideal sheaves $J$
can be computed by `strict_transform(phi, J)`.
Ideal sheaves on toric varieties are currently (30 Jan 2025) considered
experimental.


## Constructors

The main constructor returns the toric blowup induced by the star
subdivision of the polyhedral fan along a primitive vector in the
support of the fan (see 11.1 Star Subdivisions in [CLS11](@cite)):
```@docs
blow_up(X::NormalToricVarietyType, primitive_vector::AbstractVector{<:IntegerUnion}; coordinate_name::Union{String, Nothing} = nothing)
```
The ray can alternatively be given using minimal supercone coordinates:
```@docs
blow_up_along_minimal_supercone_coordinates(X::NormalToricVarietyType, minimal_supercone_coords::AbstractVector{<:RationalUnion}; coordinate_name::Union{String, Nothing} = nothing)
```
The blowup along the $n$-th (nonzero) cone in the fan of the variety $X$
is the morphism induced by the star subdivision along the barycenter of
the cone (the primitive generator of the sum of its rays):
```@docs
blow_up(X::NormalToricVarietyType, n::Int; coordinate_name::Union{String, Nothing} = nothing)
```


## Attributes

```@docs
index_of_exceptional_ray(phi::ToricBlowupMorphism)
minimal_supercone_coordinates_of_exceptional_ray(phi::ToricBlowupMorphism)
exceptional_prime_divisor(phi::ToricBlowupMorphism)
center(phi::ToricBlowupMorphism)
```
Also the following attributes of toric morphisms are supported for toric blowups:
- `lattice_homomorphism(phi::ToricBlowupMorphism)`,
- `morphism_on_torusinvariant_weil_divisor_group(phi::ToricBlowupMorphism)`,
- `morphism_on_torusinvariant_cartier_divisor_group(phi::ToricBlowupMorphism)`,
- `morphism_on_class_group(phi::ToricBlowupMorphism)`,
- `morphism_on_picard_group(phi::ToricBlowupMorphism)`.


## Methods

We can compute the total and strict transforms of homogeneous ideals in
Cox rings under star subdivisions along a primitive vector.
```@docs
strict_transform(phi::ToricBlowupMorphism, I::MPolyIdeal)
strict_transform_with_index(phi::ToricBlowupMorphism, I::MPolyIdeal)
total_transform(phi::ToricBlowupMorphism, I::Union{MPolyIdeal, MPolyDecRingElem})
```
The above functions are implemented using a $\mathbb{C}$-module
homomorphism between the Cox rings, considered as $\mathbb{C}$-modules,
that takes homogeneous ideals to homogeneous ideals:
```@docs
cox_ring_module_homomorphism
```
