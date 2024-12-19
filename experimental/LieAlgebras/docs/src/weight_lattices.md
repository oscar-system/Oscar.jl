```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Weight lattices

Weight lattices are represented by objects of type `WeightLattice <: AdditiveGroup`, and their elements by `WeightLatticeElem <: AdditiveGroupElement`.

They are introduced to have a formal parent object of all weights that correspond to a common given root system.

!!! note
    See [Cartan types](@ref) for our conventions on Cartan types and ordering of simple roots.

## Table of contents

```@contents
Pages = ["weight_lattices.md"]
Depth = 2:5
```

## Constructing weight lattices
```@docs
weight_lattice(::RootSystem)
```


## Properties of weight lattices

```@docs
rank(::WeightLattice)
is_finite(::WeightLattice)
zero(::WeightLattice)
gen(::WeightLattice, ::Int)
gens(::WeightLattice)
```

```@docs
root_system(::WeightLattice)
```


## Weight lattice elements

```@docs
WeightLatticeElem(::WeightLattice, ::Vector{<:IntegerUnion})
WeightLatticeElem(::RootSystem, ::Vector{<:IntegerUnion})
WeightLatticeElem(::WeightLattice, ::ZZMatrix)
WeightLatticeElem(::RootSystem, ::ZZMatrix)
WeightLatticeElem(::RootSpaceElem)
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
reflect(::WeightLatticeElem, ::RootSpaceElem)
reflect!(::WeightLatticeElem, ::RootSpaceElem)
```

### Conjugate dominant weight
```@docs
conjugate_dominant_weight(::WeightLatticeElem)
conjugate_dominant_weight_with_elem(::WeightLatticeElem)
```
