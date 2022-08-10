```@meta
CurrentModule = Oscar
```

```@contents
Pages = ["CohomologyClasses.md"]
```


# Cohomology Classes


## Constructors

### General constructors

```@docs
CohomologyClass(d::ToricDivisor)
CohomologyClass(c::ToricDivisorClass)
CohomologyClass(l::ToricLineBundle)
```

### Addition, subtraction and scalar multiplication

Cohomology classes can be added and subtracted via the usual `+` and `-`
operators. Moreover, multiplication by scalars from the left is supported
for scalars which are integers or of type `fmpz` or `fmpq`.

### Wedge product

The wedge product of cohomology classes is implemented via `*`.
This makes sense, since cohomology classes on (simplicial and complete)
toric varieties are elements of the cohomology ring, which then
turns out to be (a certain) quotient of the Cox ring. Hence, internally,
a cohomology class can be thought of as a polynomial in this ring and
the wedge product corresponds to the product of two (equivalence classes
of) polynomials.

A cohomology class can be wedged `n`-times with itself via `^n`,
where `n` can be an integer of of type `fmpz`.


## Properties

One can check if a cohomology class is trivial via `is_trivial`.

Equality of cohomology classes can be tested via `==`.


## Attributes

```@docs
toric_variety(c::CohomologyClass)
coefficients(c::CohomologyClass)
exponents(c::CohomologyClass)
polynomial(c::CohomologyClass)
polynomial(c::CohomologyClass, ring::MPolyQuo)
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
