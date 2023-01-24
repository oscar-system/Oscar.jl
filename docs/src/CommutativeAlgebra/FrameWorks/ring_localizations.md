```@meta
CurrentModule = Oscar 
DocTestSetup = quote
  using Oscar
end
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["localizations.md"]
```

# A Framework for Localizing Rings

For the convenience of the developer, we outline a general framework for creating concrete instances of localized rings in OSCAR,
addressing relevant abstract types as well as a standardized set of functions whose concrete behaviour must be implemented.

We roughly follow the outline of the previous subsection on localizing multivariate rings which provides illustrating examples.
With regard to notation, the name `Rloc` will refer to the localization of a commutative ring `R` with 1.

## Localized Rings

All multiplicatively closed subsets should belong to the `AbsMultSet` abstract type and all
localized rings should belong to the `AbsLocalizedRing` abstract type.

The basic functionality that has to be realized for any concrete instance of `AbsMultSet`
is the containment check for elements in multiplicatively closed subsets via the `in` function.

For each concrete instance of `AbsLocalizedRing`, the `Localization` constructor as well as the
functions `base_ring` and `inverted_set` need to be implemented. Moreover, as for any other type
of rings in OSCAR, methods for the standardized set of functions of OSCAR's
general [Ring Interface](@ref) must be supplied.

## Elements of Localized Rings

All elements of localized rings should belong to the `AbsLocalizedRingElem` abstract type.

Coercing (pairs of) elements of `R` into fractions in `Rloc` must be possible as indicated below:

```
   (Rloc::AbsLocalizedRing)(f::RingElem)
   (Rloc::AbsLocalizedRing)(f::RingElem, g::RingElem; check::Bool=true)
```
   
The first constructor maps the element `f` of `R` to the fraction `f//1` in `Rloc`.
The second constructor takes a pair `f, g` of elements of `R` to the fraction `f//g`
in `Rloc`. The default `check = true` stands for testing whether `g` is an admissible
denominator. As this test is often expensive, it may be convenient
to set `check = false`.

For any concrete instance of type `AbsLocalizedRingElem`, methods for the functions 
`parent`, `numerator`, and `denominator` must be provided. Moreover,
if a cancellation function for the type of fractions under consideration is
not yet available, such a function should be implemented and named
`reduce_fraction`.

## Homomorphisms From Localized Rings

The abstract type for homomorphisms from localized rings is `MPolyLocalizedRingHom`.
For each concrete instance, the functions `domain` and `codomain` as well as `restricted_map`
must be realized. Here, the latter function is meant to return the composition with the
localization map.

## Ideals in Localized Rings

All ideals in localized rings belong to the abstract type `AbsLocalizedIdeal`.
For a concrete instance, the constructors to be implemented are:

```
   ideal(W::AbsLocalizedRing, f::AbsLocalizedRingElem)
   ideal(W::AbsLocalizedRing, v::Vector{LocalizedRingElemType}) where {LocalizedRingElemType<:AbsLocalizedRingElem}
```

The usual getter functions  `base_ring`, `gens`, `ngens`, and `gen`   must be realized.

Moreover, a method for ideal membership via the `in` function is required.




