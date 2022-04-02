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

```@docs
Base.:+(cc1::CohomologyClass, cc2::CohomologyClass)
Base.:-(cc1::CohomologyClass, cc2::CohomologyClass)
Base.:*(c::fmpq, cc::CohomologyClass)
```

### Wedge product

```@docs
Base.:*(cc1::CohomologyClass, cc2::CohomologyClass)
Base.:^(cc::CohomologyClass, p::fmpz)
```

### Equality

```@docs
Base.:(==)(cc1::CohomologyClass, cc2::CohomologyClass)
```


## Properties

```@docs
istrivial(c::CohomologyClass)
```


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
