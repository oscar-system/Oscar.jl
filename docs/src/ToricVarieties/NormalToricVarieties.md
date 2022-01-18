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
toric_affine_space(d::Int)
del_pezzo(b::Int)
hirzebruch_surface(r::Int)
toric_projective_space(d::Int)
```

### Further Constructions

```@docs
blowup_on_ith_minimal_torus_orbit(v::AbstractNormalToricVariety, n::Int)
Base.:*(v::AbstractNormalToricVariety, w::AbstractNormalToricVariety)
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
isprojective_space(v::AbstractNormalToricVariety)
isq_gorenstein(v::AbstractNormalToricVariety)
```


## Operations for Toric Varieties

### Affine Open Covering

```@docs
affine_open_covering( v::AbstractNormalToricVariety )
```

### Characters, Weil Divisors, Cartier Divisors, Class Group and Picard Group

```@docs
cartier_divisor_group(v::AbstractNormalToricVariety)
character_lattice(v::AbstractNormalToricVariety)
class_group(v::AbstractNormalToricVariety)
map_from_cartier_divisor_group_to_torus_invariant_divisor_group(v::AbstractNormalToricVariety)
map_from_cartier_divisor_group_to_picard_group(v::AbstractNormalToricVariety)
map_from_character_to_principal_divisors(v::AbstractNormalToricVariety)
map_from_weil_divisors_to_class_group(v::AbstractNormalToricVariety)
picard_group(v::AbstractNormalToricVariety)
torusinvariant_divisor_group(v::AbstractNormalToricVariety)
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

```@docs
cox_ring(v::AbstractNormalToricVariety)
irrelevant_ideal(v::AbstractNormalToricVariety)
stanley_reisner_ideal(v::AbstractNormalToricVariety)
toric_ideal(antv::AffineNormalToricVariety)
```


## Auxillary Methods

```@docs
binomial_exponents_to_ideal(binoms::Union{AbstractMatrix, fmpz_mat})
toric_ideal(pts::Union{AbstractMatrix, fmpz_mat})
```
