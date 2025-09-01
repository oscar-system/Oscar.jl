```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```


# Cohomology Classes


## Constructors

### General constructors

Cohomology classes are elements of the cohomology ring, which is a quotient ring for normal toric varieties that
are simplicial and complete. All of the following methods accept an optional argument `completeness_check`
(default: `true`).  

By default, completeness is checked to ensure mathematical correctness. However, these checks can be computationally
expensive—especially for large toric varieties, where operations in the quotient ring (such as Gröbner basis
computations) are already costly. To speed up computations, you may disable this step by passing
`completeness_check = false`.

**Caveat:** Disabling the completeness check skips verification that the underlying toric variety is complete.
Use this option only if completeness is already known or guaranteed, and faster computations are preferred.

```@docs
cohomology_class(v::NormalToricVarietyType, p::MPolyQuoRingElem; completeness_check::Bool = true)
cohomology_class(d::ToricDivisor; completeness_check::Bool = true)
cohomology_class(c::ToricDivisorClass; completeness_check::Bool = true)
cohomology_class(l::ToricLineBundle; completeness_check::Bool = true)
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
simplify!(c::CohomologyClass)
coefficients(c::CohomologyClass)
exponents(c::CohomologyClass)
polynomial(c::CohomologyClass)
polynomial(ring::MPolyQuoRing, c::CohomologyClass)
```


## Methods

```@docs
integrate(c::CohomologyClass; completeness_check::Bool = true)
```


## Special attributes of toric varieties

```@docs
cohomology_ring(v::NormalToricVarietyType; completeness_check::Bool = true)
volume_form(v::NormalToricVariety)
intersection_form(v::NormalToricVariety)
chern_class(v::NormalToricVariety, k::Int; completeness_check::Bool = true)
chern_classes(v::NormalToricVariety; completeness_check::Bool = true)
basis_of_h4(v::NormalToricVariety; completeness_check::Bool = true)
```
