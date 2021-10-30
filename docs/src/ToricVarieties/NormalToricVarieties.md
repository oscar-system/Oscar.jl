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
CyclicQuotientSingularity
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
toric_ideal
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

Cyclic quotient singularities are quotients of $\mathbb{C}^2$ by the action of
$\mathbb{Z}/n\mathbb{Z}$ acting via 
$$\left(\begin{array}{cc}\xi & 0\\0 & \xi^q\end{array}\right)$$,
where $\xi$ is a $n$-th root of unity, and $0<q<n$ are two coprime integers.

For the notation we rely on [Chr91](@cite) and [Ste91](@cite).

!!! warning
    Note that the notion of Hirzebruch-Jung continued fraction differs from the
    commonly known continued fraction.

```@docs
continued_fraction_hirzebruch_jung(cqs::CyclicQuotientSingularity)
continued_fraction_hirzebruch_jung_to_rational(v::Vector{fmpz})
dual_continued_fraction_hirzebruch_jung(cqs::CyclicQuotientSingularity)
rational_to_continued_fraction_hirzebruch_jung(r::fmpq)
```
