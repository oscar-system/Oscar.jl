```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Weyl groups

Weyl groups are represented by objects of type `WeylGroup <: Group`, and their elements by `WeylGroupElem <: GroupElement`.

!!! warning
    Weyl groups in OSCAR afford both left and right actions on roots and weights.
    Note however, that **the left action** is the default (to align with the literature),
    and all more complex functionality is defined with respect to the left action, e.g.
    [`conjugate_dominant_weight_with_elem(::WeightLatticeElem)`](@ref).

## Table of contents

```@contents
Pages = ["weyl_groups.md"]
Depth = 2:5
```

## Constructing Weyl groups
```@docs
weyl_group(::RootSystem)
weyl_group(::ZZMatrix)
weyl_group(::Symbol, ::Int)
weyl_group(::Vector{Tuple{Symbol,Int}})
```

## Basic properties
Basic group arithmetic like `*`, and `inv` are defined for `WeylGroupElem` objects.

Using `(W::WeylGroup)(word::Vector{<:Integer})`, one can construct group elements from a word in the generators.

```@docs
is_finite(::WeylGroup)
one(::WeylGroup)
isone(::WeylGroupElem)
gen(::WeylGroup, ::Int)
gens(::WeylGroup)
number_of_generators(::WeylGroup)
order(::Type{T}, ::WeylGroup) where {T}
```

```@docs
root_system(::WeylGroup)
```

### Words and length
```@docs
word(::WeylGroupElem)
length(::WeylGroupElem)
longest_element(::WeylGroup)
```

### Bruhat order
```@docs
<(::WeylGroupElem, ::WeylGroupElem)
```


## Conversion to other group types

For many computations, it may be suitable to have a `WeylGroup` as a different kind of group object, to e.g. use functionality that is only available for that other type.

The conversion functions come in pairs: one only creates an isomorphic group object, the other also computes the isomorphism.

```@docs
fp_group(::WeylGroup)
isomorphism(::Type{FPGroup}, ::WeylGroup)
```


## Reduced expressions

```@docs
reduced_expressions(::WeylGroupElem)
```


## Action on roots and weights

```@docs
*(::WeylGroupElem, ::Union{RootSpaceElem,WeightLatticeElem})
*(::Union{RootSpaceElem,WeightLatticeElem}, ::WeylGroupElem)
```


### Orbits
TODO
