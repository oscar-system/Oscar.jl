```@meta
CurrentModule = Oscar
```

```@contents
Pages = ["Subvarieties.md"]
```


# Subvarieties

## Introduction

We focus on simplicial toric varieties. Then, any
closed subvariety is given as the vanishing set of
a homogeneous ideal in the Cox ring of the toric
variety in question (cf. proposition 5.2.4 in
[CLS11](@cite)). As of now, we provide elementary
support for closd subvarieties of simplicial toric
varieties.


## Constructors

### General constructors

```@docs
ClosedSubvarietyOfToricVariety(toric_variety::AbstractNormalToricVariety, defining_polynomials::Vector{MPolyElem_dec{fmpq, fmpq_mpoly}})
```

### Special constructors

Empty set and entire variety?


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
