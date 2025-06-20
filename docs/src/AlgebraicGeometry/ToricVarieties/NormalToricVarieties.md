```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Normal Toric Varieties

## Introduction

We introduce two main types of normal toric varieties, distinguishing between
the affine and non-affine case:
- `AffineNormalToricVariety` is the type for toric varieties associated to a cone $\sigma$, denoted by $U_{\sigma}$ in [CLS11](@cite)
- `NormalToricVariety` is the type for toric varieties associated to a polyhedral fan $\Sigma$, denoted by $X_{\Sigma}$ in [CLS11](@cite)

!!! warning
    The lattice is always assumed to be the standard lattice $\mathbb{Z}^n$.
    Transformations for non-standard lattices will have to be done by the user.


## Equality of Normal Toric Varieties

Equality `==` of normal toric varieties checks equality of the
corresponding polyhedral fans as sets of cones.
This computes the rays of both of the toric varieties, which can be
expensive if they are not already computed, meaning if
`"RAYS" in Polymake.list_properties(Oscar.pm_object(polyhedral_fan(X)))`
is false for one of the varieties.
Triple-equality `===` always checks equality of memory locations in OSCAR.


## Constructors

### Affine Toric Varieties

```@docs
affine_normal_toric_variety(C::Cone)
normal_toric_variety(C::Cone)
affine_normal_toric_variety(v::NormalToricVariety)
```

### Normal Toric Varieties

```@docs
normal_toric_variety
```

### Famous Toric Varieties

The constructors of `del_pezzo_surface`, `hirzebruch_surface`, `projective_space` and `weighted_projective_space ` *always* make a default/standard choice for the grading of the Cox ring.

```@docs
affine_space(::Type{NormalToricVariety}, d::Int)
del_pezzo_surface(::Type{NormalToricVariety}, b::Int)
hirzebruch_surface(::Type{NormalToricVariety}, r::Int)
projective_space(::Type{NormalToricVariety}, d::Int)
weighted_projective_space(::Type{NormalToricVariety}, w::Vector{T}) where {T <: IntegerUnion}
```

### Constructions based on triangulations

It is possible to associate toric varieties to star triangulations
of the lattice points of polyhedrons. Specifically, we can associate
to any full star triangulation of the lattice points of the polyhedron
in question a toric variety. For this task, we provide the following
constructors.
```@docs
normal_toric_variety_from_star_triangulation(P::Polyhedron)
normal_toric_varieties_from_star_triangulations(P::Polyhedron)
```
An application of this functionality exists in the physics.
Witten's Generalized-Sigma models (GLSM) [Wit88](@cite)
originally sparked interest in the physics community in toric varieties.
On a mathematical level, this establishes a construction of toric
varieties for which a Z^n grading of the Cox ring is provided. See
for example [FJR17](@cite), which describes this as GIT
construction [CLS11](@cite).

Explicitly, given the grading of the Cox ring, the map from
the group of torus invariant Weil divisors to the class group
is known. Under the assumption that the variety in question
has no torus factor, we can then identify the map from the
lattice to the group of torus invariant Weil divisors as the
kernel of the map from the torus invariant Weil divisor to the
class group. The latter is a map between free Abelian groups, i.e.
is provided by an integer valued matrix. The rows of this matrix
are nothing but the ray generators of the fan of the toric variety.
It then remains to triangulate these rays, hence in general for
a GLSM the toric variety is only unique up to fine regular
star triangulations. We provide the following two constructors:
```@docs
normal_toric_variety_from_glsm(charges::ZZMatrix)
normal_toric_varieties_from_glsm(charges::ZZMatrix)
```

### Further Constructions

```@docs
Base.:*(v::NormalToricVarietyType, w::NormalToricVarietyType)
projectivization(E::ToricLineBundle...)
total_space(E::ToricLineBundle...)
```


## Properties of Toric Varieties

```@docs
has_torusfactor(v::NormalToricVarietyType)
is_affine(v::NormalToricVarietyType)
is_complete(v::NormalToricVarietyType)
is_fano(v::NormalToricVarietyType)
is_gorenstein(v::NormalToricVarietyType)
is_simplicial(v::NormalToricVarietyType)
is_smooth(v::NormalToricVarietyType)
is_normal(v::NormalToricVarietyType)
is_orbifold(v::NormalToricVarietyType)
is_projective(v::NormalToricVarietyType)
is_projective_space(v::NormalToricVarietyType)
is_q_gorenstein(v::NormalToricVarietyType)
```


## Operations for Toric Varieties

### Affine Open Covering

```@docs
affine_open_covering( v::NormalToricVarietyType )
```

### Characters, Weil Divisors, Cartier Divisors, Class Group and Picard Group

```@docs
torusinvariant_cartier_divisor_group(v::NormalToricVarietyType)
character_lattice(v::NormalToricVarietyType)
class_group(v::NormalToricVarietyType)
map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(v::NormalToricVarietyType)
map_from_torusinvariant_cartier_divisor_group_to_picard_group(v::NormalToricVarietyType)
map_from_character_lattice_to_torusinvariant_weil_divisor_group(v::NormalToricVarietyType)
map_from_torusinvariant_weil_divisor_group_to_class_group(v::NormalToricVarietyType)
picard_group(v::NormalToricVarietyType)
torusinvariant_weil_divisor_group(v::NormalToricVarietyType)
torusinvariant_prime_divisors(v::NormalToricVarietyType)
```

### Gorenstein and Picard index

```@docs
gorenstein_index(v::NormalToricVarietyType)
picard_index(v::NormalToricVarietyType)
```

### Cones and Fans

```@docs
polyhedral_fan(v::NormalToricVarietyType)
cone(v::AffineNormalToricVariety)
weight_cone(v::AffineNormalToricVariety)
hilbert_basis(v::AffineNormalToricVariety)
mori_cone(v::NormalToricVariety)
nef_cone(v::NormalToricVariety)
```

### Dimensions

```@docs
dim(v::NormalToricVarietyType)
dim_of_torusfactor(v::NormalToricVarietyType)
euler_characteristic(v::NormalToricVarietyType)
betti_number(v::NormalToricVarietyType, i::Int)
```

### Rings and ideals

We support the following rings and ideals for toric varieties:
- Cox ring (also termed the "total coordinate ring" in [CLS11](@cite)),
- coordinate ring of torus,
- cohomology_ring,
- Chow ring,
- irrelevant ideal,
- Stanley-Reisner ideal,
- ideal of linear relations,
- toric ideal.
Of course, for any of these coordinate names and the coefficient ring
have to be chosen. The coefficient ring is fixed to `Q`. Therefore, the
method `coefficient_ring(v::NormalToricVarietyType)` always return
the field of rational numbers. For the coordinate names, we provide the
following setter functions:
```@docs
set_coordinate_names(v::NormalToricVarietyType, coordinate_names::AbstractVector{<:VarName})
set_coordinate_names_of_torus(v::NormalToricVarietyType, coordinate_names::AbstractVector{<:VarName})
```
The following methods allow to extract the chosen coordinates:
```@docs
coordinate_names(v::NormalToricVarietyType)
coordinate_names_of_torus(v::NormalToricVarietyType)
```
In order to efficiently construct algebraic cycles (elements of the Cox ring),
cohomology classes (elements of the cohomology ring), or in order to compare ideals,
it is imperative to fix choices of the coordinate names. The default value for
coordinate names is `[x1, x2, ... ]`. The choice of coordinate names is fixed,
once one of the above-mentioned rings is computed via one the following methods:
```@docs
cox_ring(v::NormalToricVarietyType)
irrelevant_ideal(v::NormalToricVarietyType)
ideal_of_linear_relations(v::NormalToricVarietyType)
stanley_reisner_ideal(v::NormalToricVarietyType)
toric_ideal(antv::AffineNormalToricVariety)
coordinate_ring_of_torus(v::NormalToricVarietyType)
```
One can check the status as follows:
```@docs
is_finalized(v::NormalToricVarietyType)
```
After the variety finalized, one can enforce to obtain the above ideals in different rings.
Also, one can opt to compute the above rings with a different choice of coordinate names
and different coefficient ring. To this end, one provides a custom ring (which
reflects the desired choice of coordinate names and coefficient ring) as first argument.
However, note that the cached ideals and rings are *not* altered.
```@docs
cox_ring(R::MPolyRing, v::NormalToricVarietyType)
irrelevant_ideal(R::MPolyRing, v::NormalToricVarietyType)
ideal_of_linear_relations(R::MPolyRing, v::NormalToricVarietyType)
stanley_reisner_ideal(R::MPolyRing, v::NormalToricVarietyType)
toric_ideal(R::MPolyRing, antv::AffineNormalToricVariety)
coordinate_ring_of_torus(R::MPolyRing, v::NormalToricVarietyType)
```
Along the same lines, characters can be turned into rational functions:
```@docs
character_to_rational_function(v::NormalToricVarietyType, character::Vector{ZZRingElem})
character_to_rational_function(R::MPolyRing, v::NormalToricVarietyType, character::Vector{ZZRingElem})
```


## Auxiliary Methods

```@docs
binomial_exponents_to_ideal(binoms::Union{AbstractMatrix, ZZMatrix})
toric_ideal(pts::ZZMatrix)
```
