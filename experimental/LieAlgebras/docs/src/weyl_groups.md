```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Weyl groups

```@docs
weyl_group(::RootSystem)
weyl_group(::ZZMatrix)
weyl_group(::Symbol, ::Int)
weyl_group(::Vector{Tuple{Symbol,Int}})
```

TODO: add documentation about parent call syntax `W([1,2,2])`

```@docs
is_finite(::WeylGroup)
one(::WeylGroup)
isone(::WeylGroupElem)
gen(::WeylGroup, ::Int)
gens(::WeylGroup)
number_of_generators(::WeylGroup)
order(::Type{T}, ::WeylGroup) where {T}
longest_element(::WeylGroup)
```

```@docs
root_system(::WeylGroup)
```

Basic group arithmetic like `*`, and `inv` are defined for `WeylGroupElem` objects.

```@docs
length(::WeylGroupElem)
word(::WeylGroupElem)
```

```@docs
<(::WeylGroupElem, ::WeylGroupElem)
```

```@docs
*(::WeylGroupElem, ::Union{RootSpaceElem,WeightLatticeElem})
*(::Union{RootSpaceElem,WeightLatticeElem}, ::WeylGroupElem)
```

## Conversion to other group types

For many computations, it may be suitable to have a `WeylGroup` as a different kind of group object, to e.g. use functionality that is only available for that other type.

The conversion functions come in pairs: one only creates an isomorphic group object, the other also computes the isomorphism.

```@docs
fp_group(::WeylGroup)
isomorphism(::Type{FPGroup}, ::WeylGroup)
```

## Reduced expressions
TODO
