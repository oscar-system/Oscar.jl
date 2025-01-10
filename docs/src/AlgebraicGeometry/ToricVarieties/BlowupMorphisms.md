```@meta
CurrentModule = Oscar
```

# Toric Blowups (Experimental)

It is a common goal in algebraic geometry to resolve singularities.
Certainly, toric varieties and their subvarieties are no exception and
we provide a growing set of functionality for such tasks.

In general, resolutions of toric varieties need not be toric. Indeed,
some of the functionality below requires fully-fledge schemes machinery,
which -- as of October 2023 -- is still in Oscar's experimental state.
For this reason, the methods below should be considered experimental.

We focus mainly on toric blowups given by a star subdivision of a
polyhedral fan along a primitive vector, see 11.1 Star Subdivisions in
[CLS11](@cite). Below, we refer to this new primitive vector as
`new_ray`. The main constructor is the following
- `blow_up(Y::NormalToricVariety, new_ray::AbstractVector{<:IntegerUnion}; coordinate_name::String)`
This will also construct the underlying toric morphism. We can specify
the name for the coordinate in the Cox ring that is assigned to
`new_ray` using the optional argument `coordinate_name`.

More generally, we can construct a blowup along a closed subscheme given
by an ideal in the Cox ring or by an ideal sheaf of the corresponding
covered scheme. In general, this will result in a non-toric variety.


## Constructors

The following methods blow up toric varieties. The closed subscheme
along which the blowup is constructed can be provided in different
formats. We discuss the methods in ascending generality.

For our most specialized blowup method, we focus on the n-th cone in the
fan of the variety `v` in question. This cone need not be maximal. The
ensuing star subdivision will subdivide this cone about its "diagonal"
(the sum of all its rays). The result of this will always be a toric
variety:
```@docs
blow_up(v::NormalToricVariety, n::Int; coordinate_name::String = "e")
```
More generally, we can provide a primitive element in the fan of the
variety in question and construct a toric morphism as in Section 11.1
Star Subdivisions in [CLS11](@cite). The resulting star subdivision
leads to a polyhedral fan, or put differently, the blowup is always
toric:
```@docs
blow_up(v::NormalToricVariety, new_ray::AbstractVector{<:IntegerUnion}; coordinate_name::String = "e")
```
Most generally, we encode the closed subscheme along which we blow up by
a homogeneous ideal in the Cox ring. Such blowups are often non-toric,
i.e. the return value of the following method could well be non-toric.
```@docs
blow_up(v::NormalToricVariety, I::MPolyIdeal; coordinate_name::String = "e")
```
Instead of providing the ideal, it is possible to turn a homogeneous
ideal in the Cox ring into an ideal sheaf. Consequently, we also provide
the support for the following method.
```@docs
blow_up(m::NormalToricVariety, I::ToricIdealSheafFromCoxRingIdeal; coordinate_name::String = "e")
```


## Attributes

```@docs
underlying_morphism(bl::ToricBlowupMorphism)
index_of_new_ray(bl::ToricBlowupMorphism)
center_data(bl::ToricBlowupMorphism)
center_unnormalized(bl::ToricBlowupMorphism)
exceptional_prime_divisor(bl::ToricBlowupMorphism)
```
Based on `underlying_morphism`, also the following attributes of toric
morphisms are supported for toric blowups:
- `grid_morphism(bl::ToricBlowupMorphism)`,
- `morphism_on_torusinvariant_weil_divisor_group(bl::ToricBlowupMorphism)`,
- `morphism_on_torusinvariant_cartier_divisor_group(bl::ToricBlowupMorphism)`,
- `morphism_on_class_group(bl::ToricBlowupMorphism)`,
- `morphism_on_picard_group(bl::ToricBlowupMorphism)`.
The total and strict transform of ideal sheaves along blowups, not
necessarily toric, can be computed:
```@docs
total_transform(f::AbsSimpleBlowupMorphism, II::AbsIdealSheaf)
```


## Arithmetics

Toric blowups can be added, subtracted and multiplied by rational
numbers. The results of such operations will not be toric morphisms,
i.e. they no longer correspond to the blowup of a certain closed
subscheme. Arithmetics among toric blowups and general toric morphisms
is also supported, as well as equality for toric blowups.
