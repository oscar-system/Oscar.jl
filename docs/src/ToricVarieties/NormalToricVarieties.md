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
hastorusfactor(v::AbstractNormalToricVariety)
isaffine(v::AbstractNormalToricVariety)
iscomplete(v::AbstractNormalToricVariety)
isfano(v::AbstractNormalToricVariety)
isgorenstein(v::AbstractNormalToricVariety)
issimplicial(v::AbstractNormalToricVariety)
issmooth(v::AbstractNormalToricVariety)
isnormal(v::AbstractNormalToricVariety)
isorbifold(v::AbstractNormalToricVariety)
isprojective(v::AbstractNormalToricVariety)
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

### Coefficient ring and coordinate names

We support the Cox ring (also termed the "total coordinate ring" in [CLS11](@cite)) and the following ideals:
- `irrelevant ideal`,
- `Stanley-Reisner ideal`
For their computation, names for the indeterminates of the Cox ring must be chosen.
The default value is `[x1, x2, ... ]`. the user can overwrite this at any time.

```@docs
coordinate_names(v::AbstractNormalToricVariety)
set_coordinate_names(v::AbstractNormalToricVariety, coordinate_names::Vector{String})
```

Likewise, a coefficient ring for the Cox ring must be chosen. The default value is 
the field of rational numbers. The user can overwrite this at any time.

```@docs
coefficient_ring(v::AbstractNormalToricVariety)
set_coefficient_ring(v::AbstractNormalToricVariety, coefficient_ring::AbstractAlgebra.Ring)
```

Similarly, the computation of the coordinate ring of the torus requires coordinate names.
Again, the default value is `[x1, x2, ... ]` and be overwritten at any time.

```@docs
coordinate_names_of_torus(v::AbstractNormalToricVariety)
set_coordinate_names_of_torus(v::AbstractNormalToricVariety, coordinate_names::Vector{String})
```


### Rings and ideals

```@docs
cox_ring(v::AbstractNormalToricVariety)
cox_ring(R::MPolyRing, v::AbstractNormalToricVariety)
irrelevant_ideal(v::AbstractNormalToricVariety)
irrelevant_ideal(R::MPolyRing, v::AbstractNormalToricVariety)
stanley_reisner_ideal(v::AbstractNormalToricVariety)
stanley_reisner_ideal(R::MPolyRing, v::AbstractNormalToricVariety)
toric_ideal(antv::AffineNormalToricVariety)
toric_ideal(R::MPolyRing, antv::AffineNormalToricVariety)
coordinate_ring_of_torus(v::AbstractNormalToricVariety)
coordinate_ring_of_torus(R::MPolyRing, v::AbstractNormalToricVariety)
character_to_rational_function(v::AbstractNormalToricVariety, character::Vector{fmpz})
character_to_rational_function(R::MPolyRing, v::AbstractNormalToricVariety, character::Vector{fmpz})
```

### Sheaves

```@docs
StructureSheaf(v::AbstractNormalToricVariety)
```

## Auxillary Methods

```@docs
binomial_exponents_to_ideal(binoms::Union{AbstractMatrix, fmpz_mat})
toric_ideal(pts::fmpz_mat)
```
