```@meta
CurrentModule = Oscar
```

# Toric Blowdown Morphisms (Experimental)

It is a common goal in algebraic geometry to resolve singularities. Certainly, (sub)varieties of
toric varieties are no exception and we provide a growing set of functionality for such tasks.

In general, resolutions need not be toric. Indeed, some of the functionality below requires
fully-fledge schemes machinery, which -- as of this writing (October 2023) -- is still in
Oscar's experimental state. For this reason, the methods below should be considered experimental.


## Constructors

Blowups of toric varieties are obtained from star subdivisions of polyhedral fans. In the most general form,
a star subdivision is defined by a new primitive element in the fan. Below, we refer to this new primitive
element as `new_ray`. In addition to this `new_ray`, our design of toric blowdown morphisms requires an
underlying toric morphism. With an eye towards covered schemes as possible return value, any toric blowdown
morphism must also know (to compute) its blowup center in the form of an ideal sheaf. The following constructor
allows to set this ideal sheaf center upon construction:
- `toric_blowdown_morphism(bl::ToricMorphism, new_ray::AbstractVector{<:IntegerUnion}, center::IdealSheaf)`
The "working-horse" constructor however is the following:
- `toric_blowdown_morphism(Y::NormalToricVariety, new_ray::AbstractVector{<:IntegerUnion}, coordinate_name::String, set_attributes::Bool)`
This constructor will, among others, construct the underlying toric morphism. In addition, we can then specify
a name for the coordinate in the Cox ring that is assigned to `new_ray`.



## Blowdown morphisms from blowing up toric varieties

The following methods blow up toric varieties. The center of the blowup can be provided in different formats.
We discuss the methods in ascending generality.

For our most specialized blow-up method, we focus on the n-th cone in the fan of the variety `v` in question.
This cone need not be maximal! The ensuing star subdivision will subdivide this cone about its "diagonal"
(the sum of all its rays). The result of this will always be a toric variety:
```@docs
blow_up(v::NormalToricVarietyType, n::Int; coordinate_name::String = "e", set_attributes::Bool = true)
```
More generally, we can provide a primitive element in the fan of the variety in question. The resulting
star subdivision leads to a polyhedral fan, or put differently, the blow-up along this center is always toric:
```@docs
blow_up(v::NormalToricVarietyType, new_ray::AbstractVector{<:IntegerUnion}; coordinate_name::String = "e", set_attributes::Bool = true)
```
Finally, and most generally, we encode the blowup center by a homogeneous ideal in the Cox ring.
Such blowups center may easily lead to non-toric blowups, i.e. the return value of the following method
could well be non-toric.
```@docs
blow_up(v::NormalToricVarietyType, I::MPolyIdeal; coordinate_name::String = "e", set_attributes::Bool = true)
```


## Attributes

```@docs
underlying_morphism(bl::ToricBlowdownMorphism)
index_of_new_ray(bl::ToricBlowdownMorphism)
center(bl::ToricBlowdownMorphism)
exceptional_divisor(bl::ToricBlowdownMorphism)
```
Based on `underlying_morphism`, also the following attributes of toric morphisms are supported for toric
blowdown morphisms:
- `grid_morphism(bl::ToricBlowdownMorphism)`,
- `morphism_on_torusinvariant_weil_divisor_group(bl::ToricBlowdownMorphism)`,
- `morphism_on_torusinvariant_cartier_divisor_group(bl::ToricBlowdownMorphism)`,
- `morphism_on_class_group(bl::ToricBlowdownMorphism)`,
- `morphism_on_picard_group(bl::ToricBlowdownMorphism)`.
The total and strict transform of ideal sheaves along toric blowdown morphisms can be computed:
```@docs
total_transform(f::AbsSimpleBlowdownMorphism, II::IdealSheaf)
```


## Arithmetics

Toric blowdown morphisms can be added, subtracted and multiplied by rational numbers. The results of such
operations will be toric morphisms, i.e. no longer attributed to the blowup of a certain locus. Arithmetics
among toric blowdown morphisms and general toric morphisms is also supported, as well as equality for toric
blowdown morphisms.
