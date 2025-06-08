# Weyl groups

Weyl groups are represented by objects of type `WeylGroup <: Group`, and their elements by `WeylGroupElem <: GroupElement`.

!!! warning
    Weyl groups in OSCAR only afford **right** actions on roots and weights.
    Note however, that this may differ from the literature, but is to stay
    consistent with the conventions in the rest of OSCAR.

!!! note
    See [Cartan types](@ref) for our conventions on Cartan types and ordering of simple roots.

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

Finite Weyl groups support iteration over all group elements (in an arbitrary order).

```@docs
is_finite(::WeylGroup)
one(::WeylGroup)
isone(::WeylGroupElem)
gen(::WeylGroup, ::Int)
gens(::WeylGroup)
number_of_generators(::WeylGroup)
order(::Type{T}, ::WeylGroup) where {T}
is_finite_order(::WeylGroupElem)
order(::Type{T}, ::WeylGroupElem) where {T}
```

```@docs
root_system(::WeylGroup)
```

### Element constructors

Using `(W::WeylGroup)(word::Vector{<:Integer})`, one can construct group elements from a word in the generators.

```@docs; canonical=false
gen(::WeylGroup, ::Int)
gens(::WeylGroup)
```

```@docs
reflection(::RootSpaceElem)
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


## Reduced expressions

```@docs
reduced_expressions(::WeylGroupElem)
```


## Action on roots and weights

```@docs
*(::Union{RootSpaceElem,WeightLatticeElem}, ::WeylGroupElem)
```

```@docs
geometric_representation(::WeylGroup)
dual_geometric_representation(::WeylGroup)
```


### Orbits

```@docs
weyl_orbit(::WeightLatticeElem)
```
