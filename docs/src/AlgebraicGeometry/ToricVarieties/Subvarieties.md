```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```


# Subvarieties

## Introduction

We focus on simplicial toric varieties. Then, any
closed subvariety is given as the vanishing set of
a homogeneous ideal in the Cox ring of the toric
variety in question (cf. proposition 5.2.4 in
[CLS11](@cite)). As of now, we provide elementary
support for closed subvarieties of simplicial toric
varieties.


## Constructors

### General constructors

```@docs
closed_subvariety_of_toric_variety(toric_variety::NormalToricVarietyType, defining_polynomials::Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}})
closed_subvariety_of_toric_variety(toric_variety::NormalToricVarietyType, defining_ideal::MPolyIdeal)
```


## Properties

```@docs
is_empty(c::ClosedSubvarietyOfToricVariety)
```


## Attributes

```@docs
toric_variety(c::ClosedSubvarietyOfToricVariety)
defining_ideal(c::ClosedSubvarietyOfToricVariety)
radical(c::ClosedSubvarietyOfToricVariety)
```
