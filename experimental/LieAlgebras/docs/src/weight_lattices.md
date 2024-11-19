```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Weight lattices

TODO


## Table of contents

```@contents
Pages = ["weight_lattices.md"]
Depth = 2:5
```

## Constructing weight lattices

TODO



## Properties of weight lattices

TODO




## Weight lattice elements

```@docs
WeightLatticeElem(::RootSystem, ::Vector{<:IntegerUnion})
WeightLatticeElem(::RootSystem, ::ZZMatrix)
WeightLatticeElem(::RootSpaceElem)
zero(::Type{WeightLatticeElem}, ::RootSystem)
```

```@docs
root_system(::WeightLatticeElem)
```

Basic arithmetic operations like `zero`, `+`, `-`, `*` (with integer scalars), and `==` are supported.

```@docs
coeff(::WeightLatticeElem, ::Int)
coefficients(::WeightLatticeElem)
```

```@docs
iszero(::WeightLatticeElem)
is_dominant(::WeightLatticeElem)
is_fundamental_weight(::WeightLatticeElem)
is_fundamental_weight_with_index(::WeightLatticeElem)
```

### Reflections
```@docs
reflect(::WeightLatticeElem, ::Int)
reflect!(::WeightLatticeElem, ::Int)
```

### Conjugate dominant weight
```@docs
conjugate_dominant_weight(::WeightLatticeElem)
conjugate_dominant_weight_with_left_elem(::WeightLatticeElem)
conjugate_dominant_weight_with_right_elem(::WeightLatticeElem)
```
