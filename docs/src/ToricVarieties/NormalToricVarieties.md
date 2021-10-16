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

```@docs
AffineNormalToricVariety(C::Cone)
CyclicQuotientSingularity(n::Int64, q::Int64)
NormalToricVariety(C::Cone)
NormalToricVariety(PF::PolyhedralFan)
NormalToricVariety(P::Polyhedron)
toric_projective_space
hirzebruch_surface
del_pezzo
```


## Properties of toric varieties

```@docs
hastorusfactor
isaffine
iscomplete
isfano
isgorenstein
issimplicial
issmooth
isnormal
isorbifold
isprojective
isq_gorenstein
```


## Operations for toric varieties

### Dimensions

```@docs
dim
dim_of_torusfactor
ith_betti_number(v::AbstractNormalToricVariety, i::Int)
euler_characteristic
```

### Rings and ideals

```@docs
cox_ring
stanley_reisner_ideal
irrelevant_ideal
```

### Characters, Weil divisor and the class group

```@docs
character_lattice
torusinvariant_divisor_group
class_group
map_from_character_to_principal_divisors
map_from_weil_divisors_to_class_group
torusinvariant_prime_divisors
```

### Cones and fans

```@docs
fan_of_variety
fan
nef_cone
mori_cone
```

### Affine covering

```@docs
affine_open_covering( v::AbstractNormalToricVariety )
```

### Toric ideal

To come very soon.

### Advanced constructions

```@docs
blowup_on_ith_minimal_torus_orbit(v::AbstractNormalToricVariety, n::Int)
Base.:*(v::AbstractNormalToricVariety, w::AbstractNormalToricVariety)
```

### Comparison

```@docs
isprojective_space(v::AbstractNormalToricVariety)
```


## Cyclic Quotient Singularities
```@docs
continued_fraction(cqs::CyclicQuotientSingularity)
dual_continued_fraction(cqs::CyclicQuotientSingularity)
```
