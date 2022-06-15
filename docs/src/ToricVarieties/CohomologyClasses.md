```@meta
CurrentModule = Oscar
```

```@contents
Pages = ["CohomologyClasses.md"]
```


# Cohomology Classes


## Constructors

### Generic constructors

```@docs
CohomologyClass(d::ToricDivisor)
CohomologyClass(c::ToricDivisorClass)
CohomologyClass(l::ToricLineBundle)
```

### Addition, subtraction and scalar multiplication

Addition of cohomology classes `cc1` and `cc2` is implemented by
`cc1+cc2`. Similarly, we can subtract the classes by `cc1-cc2`.
Scalar multiplication with `c` (this could be an integer,
fmpz or even fmpq number) is supported by `c*cc1`.

### Wedge product

The wedge product of cohomology classes `cc1` and `cc2`
is computed by `cc1*cc2`. This makes sense, since cohomology
classes on toric varieties are elements of the cohomology ring, which
in turn is (a certain) quotient of the Cox ring. Hence, internally,
a cohomology class is just a polynomial in this ring and the wedge
product corresponds to the product of two (equivalence classes of)
polynomials. We also support `cc1^n`, which corresponds to
computing the wedge product of `cc1` with itself `n`-times.

### Equality

Equality of cohomology classes `cc1` and `cc2` is
implemented by `cc1 == cc2`.


## Properties

To check if a cohomology class `c` is trivial, one can invoke `is_trivial(c)`.


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
