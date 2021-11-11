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

### Affine toric varieties

```@docs
AffineNormalToricVariety(C::Cone)
NormalToricVariety(C::Cone)
AffineNormalToricVariety(v::NormalToricVariety)
```

### Cyclic Quotient Singularities

Cyclic quotient singularities are quotients of $\mathbb{C}^2$ by the action of
$\mathbb{Z}/n\mathbb{Z}$ acting via 
$$\left(\begin{array}{cc}\xi & 0\\0 & \xi^q\end{array}\right)$$,
where $\xi$ is a $n$-th root of unity, and $0<q<n$ are two coprime integers.

For the notation we rely on [Chr91](@cite) and [Ste91](@cite).

!!! warning
    Note that the notion of Hirzebruch-Jung continued fraction differs from the
    commonly known continued fraction.

```@docs
CyclicQuotientSingularity(n::fmpz, q::fmpz)
continued_fraction_hirzebruch_jung(cqs::CyclicQuotientSingularity)
continued_fraction_hirzebruch_jung_to_rational(v::Vector{fmpz})
dual_continued_fraction_hirzebruch_jung(cqs::CyclicQuotientSingularity)
rational_to_continued_fraction_hirzebruch_jung(r::fmpq)
```

### Normal toric varieties

```@docs
NormalToricVariety(PF::PolyhedralFan)
NormalToricVariety(P::Polyhedron)
```

### Famous toric vareties

```@docs
del_pezzo(b::Int)
hirzebruch_surface(r::Int)
toric_projective_space(d::Int)
```

### Further constructions

```@docs
blowup_on_ith_minimal_torus_orbit(v::AbstractNormalToricVariety, n::Int)
Base.:*(v::AbstractNormalToricVariety, w::AbstractNormalToricVariety)
```


## Properties of toric varieties

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


## Operations for toric varieties

### Affine open covering

```@docs
affine_open_covering( v::AbstractNormalToricVariety )
```

### Characters, Weil divisors and the class group

```@docs
character_lattice(v::AbstractNormalToricVariety)
class_group(v::AbstractNormalToricVariety)
map_from_character_to_principal_divisors(v::AbstractNormalToricVariety)
map_from_weil_divisors_to_class_group(v::AbstractNormalToricVariety)
torusinvariant_divisor_group(v::AbstractNormalToricVariety)
torusinvariant_prime_divisors(v::AbstractNormalToricVariety)
```

### Cones and fans

```@docs
fan_of_variety(v::NormalToricVariety)
fan(v::NormalToricVariety)
mori_cone(v::NormalToricVariety)
nef_cone(v::NormalToricVariety)
```

### Dimensions

```@docs
dim(v::AbstractNormalToricVariety)
dim_of_torusfactor(v::AbstractNormalToricVariety)
euler_characteristic(v::AbstractNormalToricVariety)
ith_betti_number(v::AbstractNormalToricVariety, i::Int)
```

### Rings and ideals

```@docs
cox_ring(v::AbstractNormalToricVariety)
irrelevant_ideal(v::AbstractNormalToricVariety)
stanley_reisner_ideal(v::AbstractNormalToricVariety)
```

### Toric ideals

```@docs
binomial_exponents_to_ideal(binoms::Union{AbstractMatrix, fmpz_mat})
toric_ideal(pts::Union{AbstractMatrix, fmpz_mat})
toric_ideal(antv::AffineNormalToricVariety)
```
