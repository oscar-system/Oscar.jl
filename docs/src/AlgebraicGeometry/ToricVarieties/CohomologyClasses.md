```@meta
CurrentModule = Oscar
```


# Cohomology Classes


## Constructors

### General constructors

```@docs
cohomology_class(v::AbstractNormalToricVariety, p::MPolyQuoRingElem)
cohomology_class(d::ToricDivisor)
cohomology_class(c::ToricDivisorClass)
cohomology_class(l::ToricLineBundle)
```

### Addition, subtraction and scalar multiplication

Cohomology classes can be added and subtracted via the usual `+` and `-`
operators. Moreover, multiplication by scalars from the left is supported
for scalars which are integers or of type `ZZRingElem` or `QQFieldElem`.

### Wedge product

The wedge product of cohomology classes is implemented via `*`, 
using internally the multiplication of the corresponding polynomial 
(equivalence classes) in the Cox ring.

A cohomology class can be wedged `n`-times with itself via `^n`,
where `n` can be an integer or of type `ZZRingElem`.


## Properties

One can check if a cohomology class is trivial via `is_trivial`.

Equality of cohomology classes can be tested via `==`.


## Attributes

```@docs
toric_variety(c::CohomologyClass)
coefficients(c::CohomologyClass)
exponents(c::CohomologyClass)
polynomial(c::CohomologyClass)
polynomial(ring::MPolyQuoRing, c::CohomologyClass)
```


## Methods

```@docs
integrate(c::CohomologyClass)
```


## Special attributes of toric varieties

```@docs
cohomology_ring(v::AbstractNormalToricVariety)
volume_form(v::NormalToricVariety)
intersection_form(v::NormalToricVariety)
```
