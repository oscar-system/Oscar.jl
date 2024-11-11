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
weyl_group(::Tuple{Symbol,Int}...)
```

TODO: add documentation about parent call syntax `W([1,2,2])`

```@docs
is_finite(::WeylGroup)
one(::WeylGroup)
gen(::WeylGroup, ::Int)
gens(::WeylGroup)
number_of_generators(::WeylGroup)
order(::Type{T}, ::WeylGroup) where {T}
longest_element(::WeylGroup)
```

```@docs
root_system(::WeylGroup)
```

```@docs
*(::WeylGroupElem, ::WeylGroupElem)
inv(::WeylGroupElem)
isone(::WeylGroupElem)
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
