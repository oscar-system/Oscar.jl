```@meta
CurrentModule = Oscar
```

```@contents
Pages = ["NormalToricVarieties.md"]
```

# Normal Toric Varieties

## Introduction

We introduce two main types of normal toric varieties, distinguishing between
the affine and non-affine case:
- `AffineNormalToricVariety` is the toric variety associated to a cone $\sigma$, denoted by $U_{\sigma}$ in [CLS11](@cite)
- `NormalToricVariety` is the toric variety associated to a polyhedral fan $\Sigma$, denoted by $X_{\Sigma}$ in [CLS11](@cite)

!!! warning
    The lattice is always assumed to be the standard lattice $\mathbb{Z}^n$.
    Transformations for non-standard lattices will have to be done by the user.


## Constructors

### Affine Toric Varieties

```@docs
AffineNormalToricVariety(C::Cone)
NormalToricVariety(C::Cone)
AffineNormalToricVariety(v::NormalToricVariety)
```

### Normal Toric Varieties

```@docs
NormalToricVariety(rays::Vector{Vector{Int64}}, max_cones::Vector{Vector{Int64}})
NormalToricVariety(PF::PolyhedralFan)
NormalToricVariety(P::Polyhedron)
```

### Famous Toric Vareties

```@docs
affine_space(::Type{NormalToricVariety}, d::Int)
del_pezzo(b::Int)
hirzebruch_surface(r::Int)
projective_space(::Type{NormalToricVariety}, d::Int)
```

### Further Constructions

```@docs
blowup_on_ith_minimal_torus_orbit(v::AbstractNormalToricVariety, n::Int, coordinate_name::String)
Base.:*(v::AbstractNormalToricVariety, w::AbstractNormalToricVariety)
NormalToricVarietiesFromStarTriangulations(P::Polyhedron)
NormalToricVarietyFromGLSM(charges::fmpz_mat)
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

### Cones and Fans

```@docs
fan(v::AbstractNormalToricVariety)
cone(v::AffineNormalToricVariety)
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
character_to_rational_function(v::AbstractNormalToricVariety, character::Vector{fmpz})
character_to_rational_function(R::MPolyRing, v::AbstractNormalToricVariety, character::Vector{fmpz})
```


## Auxillary Methods

```@docs
binomial_exponents_to_ideal(binoms::Union{AbstractMatrix, fmpz_mat})
toric_ideal(pts::fmpz_mat)
```
