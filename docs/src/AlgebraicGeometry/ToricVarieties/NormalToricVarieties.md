```@meta
CurrentModule = Oscar
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


## Constructors

All of the constructors below accept as their last argument an optional boolean. This boolean `set_attributes` controls whether the constructors
sets attributes for the constructed variety. The benefit of such setters is increased performance. However, as per usual, a shortcut comes at
a price. In the case at hand, it might lead to possible inconsistencies, despite our best efforts to prevent this from happening.

The default is `set_attributes = true`, that is our constructors set attributes upon construction. If not desired, it can be switched off by
passing `set_attributes = false` as last argument. Note however, that the constructors of `del_pezzo_surface`, `hirzebruch_surface`,
`projective_space` and `weighted_projective_space ` *always* make a default/standard choice for the grading of the Cox ring.


### Affine Toric Varieties

```@docs
affine_normal_toric_variety(C::Cone; set_attributes::Bool = true)
normal_toric_variety(C::Cone; set_attributes::Bool = true)
affine_normal_toric_variety(v::NormalToricVariety; set_attributes::Bool = true)
```

### Normal Toric Varieties

```@docs
normal_toric_variety(rays::AbstractCollection[RayVector], max_cones::Vector{Vector{Int64}}; non_redundant::Bool = false)
normal_toric_variety(PF::PolyhedralFan)
normal_toric_variety(P::Polyhedron; set_attributes::Bool = true)
```

### Famous Toric Vareties

```@docs
affine_space(::Type{NormalToricVariety}, d::Int; set_attributes::Bool = true)
del_pezzo_surface(::Type{NormalToricVariety}, b::Int; set_attributes::Bool = true)
hirzebruch_surface(::Type{NormalToricVariety}, r::Int; set_attributes::Bool = true)
projective_space(::Type{NormalToricVariety}, d::Int; set_attributes::Bool = true)
weighted_projective_space(::Type{NormalToricVariety}, w::Vector{T}; set_attributes::Bool = true) where {T <: IntegerUnion}
```

### Constructions based on triangulations

It is possible to associate toric varieties to star triangulations
of the lattice points of polyhedrons. Specifically, we can associate
to any full star triangulation of the lattice points of the polyhedron
in question a toric variety. For this task, we provide the following
constructors.
```@docs
normal_toric_variety_from_star_triangulation(P::Polyhedron; set_attributes::Bool = true)
normal_toric_varieties_from_star_triangulations(P::Polyhedron; set_attributes::Bool = true)
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
normal_toric_variety_from_glsm(charges::ZZMatrix; set_attributes::Bool = true)
normal_toric_varieties_from_glsm(charges::ZZMatrix; set_attributes::Bool = true)
```

### Further Constructions

```@docs
blow_up(v::AbstractNormalToricVariety, I::MPolyIdeal; coordinate_name::String = "e", set_attributes::Bool = true)
blow_up(v::AbstractNormalToricVariety, new_ray::AbstractVector{<:IntegerUnion}; coordinate_name::String = "e", set_attributes::Bool = true)
blow_up(v::AbstractNormalToricVariety, n::Int; coordinate_name::String = "e", set_attributes::Bool = true)
Base.:*(v::AbstractNormalToricVariety, w::AbstractNormalToricVariety; set_attributes::Bool = true)
```


## Properties of Toric Varieties

```@docs
has_torusfactor(v::AbstractNormalToricVariety)
is_affine(v::AbstractNormalToricVariety)
is_complete(v::AbstractNormalToricVariety)
is_fano(v::AbstractNormalToricVariety)
is_gorenstein(v::AbstractNormalToricVariety)
is_simplicial(v::AbstractNormalToricVariety)
is_smooth(v::AbstractNormalToricVariety)
is_normal(v::AbstractNormalToricVariety)
is_orbifold(v::AbstractNormalToricVariety)
is_projective(v::AbstractNormalToricVariety)
is_projective_space(v::AbstractNormalToricVariety)
is_q_gorenstein(v::AbstractNormalToricVariety)
```


## Operations for Toric Varieties

### Affine Open Covering

```@docs
affine_open_covering( v::AbstractNormalToricVariety )
```

### Characters, Weil Divisors, Cartier Divisors, Class Group and Picard Group

```@docs
torusinvariant_cartier_divisor_group(v::AbstractNormalToricVariety)
character_lattice(v::AbstractNormalToricVariety)
class_group(v::AbstractNormalToricVariety)
map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(v::AbstractNormalToricVariety)
map_from_torusinvariant_cartier_divisor_group_to_picard_group(v::AbstractNormalToricVariety)
map_from_character_lattice_to_torusinvariant_weil_divisor_group(v::AbstractNormalToricVariety)
map_from_torusinvariant_weil_divisor_group_to_class_group(v::AbstractNormalToricVariety)
picard_group(v::AbstractNormalToricVariety)
torusinvariant_weil_divisor_group(v::AbstractNormalToricVariety)
torusinvariant_prime_divisors(v::AbstractNormalToricVariety)
```

### Gorenstein and Picard index

```@docs
gorenstein_index(v::AbstractNormalToricVariety)
picard_index(v::AbstractNormalToricVariety)
```

### Cones and Fans

```@docs
polyhedral_fan(v::AbstractNormalToricVariety)
cone(v::AffineNormalToricVariety)
dual_cone(v::AffineNormalToricVariety)
hilbert_basis(v::AffineNormalToricVariety)
mori_cone(v::NormalToricVariety)
nef_cone(v::NormalToricVariety)
```

### Dimensions

```@docs
dim(v::AbstractNormalToricVariety)
dim_of_torusfactor(v::AbstractNormalToricVariety)
euler_characteristic(v::AbstractNormalToricVariety)
betti_number(v::AbstractNormalToricVariety, i::Int)
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
method `coefficient_ring(v::AbstractNormalToricVariety)` always return
the field of rational numbers. For the coordinate names, we provide the
following setter functions:
```@docs
set_coordinate_names(v::AbstractNormalToricVariety, coordinate_names::Vector{String})
set_coordinate_names_of_torus(v::AbstractNormalToricVariety, coordinate_names::Vector{String})
```
The following methods allow to etract the chosen coordinates:
```@docs
coordinate_names(v::AbstractNormalToricVariety)
coordinate_names_of_torus(v::AbstractNormalToricVariety)
```
In order to efficiently construct algebraic cycles (elements of the Chox ring),
cohomology classes (elements of the cohomology ring), or in order to compare ideals,
it is imperative to fix choices of the coordinate names. The default value for
coordinate names is `[x1, x2, ... ]`. The choice of coordinate names is fixed,
once one of the above-mentioned rings is computed via one the following methods:
```@docs
cox_ring(v::AbstractNormalToricVariety)
irrelevant_ideal(v::AbstractNormalToricVariety)
ideal_of_linear_relations(v::AbstractNormalToricVariety)
stanley_reisner_ideal(v::AbstractNormalToricVariety)
toric_ideal(antv::AffineNormalToricVariety)
coordinate_ring_of_torus(v::AbstractNormalToricVariety)
```
One can check the status as follows:
```@docs
is_finalized(v::AbstractNormalToricVariety)
```
After the variety finalized, one can enforce to obtain the above ideals in different rings.
Also, one can opt to compute the above rings with a different choice of coordinate names
and different coefficient ring. To this end, onc provides a custom ring (which
reflects the desired choice of coordinate names and coefficient ring) as first argument.
However, note that the cached ideals and rings are *not* altered.
```@docs
cox_ring(R::MPolyRing, v::AbstractNormalToricVariety)
irrelevant_ideal(R::MPolyRing, v::AbstractNormalToricVariety)
ideal_of_linear_relations(R::MPolyRing, v::AbstractNormalToricVariety)
stanley_reisner_ideal(R::MPolyRing, v::AbstractNormalToricVariety)
toric_ideal(R::MPolyRing, antv::AffineNormalToricVariety)
coordinate_ring_of_torus(R::MPolyRing, v::AbstractNormalToricVariety)
```
Along the same lines, characters can be turned into rational functions:
```@docs
character_to_rational_function(v::AbstractNormalToricVariety, character::Vector{ZZRingElem})
character_to_rational_function(R::MPolyRing, v::AbstractNormalToricVariety, character::Vector{ZZRingElem})
```


## Auxiliary Methods

```@docs
binomial_exponents_to_ideal(binoms::Union{AbstractMatrix, ZZMatrix})
toric_ideal(pts::ZZMatrix)
```
