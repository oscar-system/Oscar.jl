```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Root systems

Root systems in this module are meant to be abstract root systems, i.e. they are represented by a set of roots (vectors in an euclidean space).

The relevant types around root systems are:
- `RootSystem` for the root system itself,
- `RootSpaceElem` for elements in the root space, i.e. roots and linear combinations thereof,
- `DualRootSpaceElem` for elements in the dual root space, i.e. coroots and linear combinations thereof.

!!! warning
    Most functionality around root systems is currently only intended to be used with root systems of finite type.
    For root systems of affine type, some documentation may be ill-phrased or incorrect, and some functions may not work as intended.

## Table of contents

```@contents
Pages = ["root_systems.md"]
Depth = 2:5
```

## Constructing root systems

```@docs
root_system(::ZZMatrix)
root_system(::Symbol, ::Int64)
root_system(::Vector{Tuple{Symbol, Int64}})
```


## Properties of root systems

```@docs
is_simple(::RootSystem)
rank(::RootSystem)
```


### Cartan matrix and Weyl group
```@docs
cartan_matrix(::RootSystem)
```
```@docs; canonical=false
weight_lattice(::RootSystem)
weyl_group(::RootSystem)
```


### Root system type
```@docs
has_root_system_type(::RootSystem)
root_system_type(::RootSystem)
root_system_type_with_ordering(::RootSystem)
```


### Root getters
```@docs
number_of_roots(::RootSystem)
number_of_positive_roots(::RootSystem)
number_of_simple_roots(::RootSystem)
```

The following functions return roots, see [Root space elements](@ref) for more information.
```@docs
root(::RootSystem, ::Int64)
roots(::RootSystem)
simple_root(::RootSystem, ::Int64)
simple_roots(::RootSystem)
positive_root(::RootSystem, ::Int64)
positive_roots(::RootSystem)
negative_root(::RootSystem, ::Int64)
negative_roots(::RootSystem)
```


### Coroot getters
The following functions return coroots, see [Dual root space elements](@ref) for more information.
```@docs
coroot(::RootSystem, ::Int64)
coroots(::RootSystem)
simple_coroot(::RootSystem, ::Int64)
simple_coroots(::RootSystem)
positive_coroot(::RootSystem, ::Int64)
positive_coroots(::RootSystem)
negative_coroot(::RootSystem, ::Int64)
negative_coroots(::RootSystem)
```


### Weight getters
The following functions return weights, see [Weight lattice elements](@ref) for more information.
```@docs
fundamental_weight(::RootSystem, ::Int64)
fundamental_weights(::RootSystem)
weyl_vector(::RootSystem)
```


## Root space elements

```@docs
RootSpaceElem(::RootSystem, ::Vector{<:RationalUnion})
RootSpaceElem(::RootSystem, ::QQMatrix)
RootSpaceElem(::WeightLatticeElem)
zero(::Type{RootSpaceElem}, ::RootSystem)
```

```@docs
root_system(::RootSpaceElem)
```

Basic arithmetic operations like `zero`, `+`, `-`, `*` (with rational scalars), and `==` are supported.

```@docs
coeff(::RootSpaceElem, ::Int)
coefficients(::RootSpaceElem)
```

```@docs
height(::RootSpaceElem)
iszero(::RootSpaceElem)
```

### Root testing
```@docs
is_root(::RootSpaceElem)
is_root_with_index(::RootSpaceElem)
is_simple_root(::RootSpaceElem)
is_simple_root_with_index(::RootSpaceElem)
is_positive_root(::RootSpaceElem)
is_positive_root_with_index(::RootSpaceElem)
is_negative_root(::RootSpaceElem)
is_negative_root_with_index(::RootSpaceElem)
```

### Reflections
```@docs
reflect(::RootSpaceElem, ::Int)
reflect!(::RootSpaceElem, ::Int)
```


## Dual root space elements

```@docs
DualRootSpaceElem(::RootSystem, ::Vector{<:RationalUnion})
DualRootSpaceElem(::RootSystem, ::QQMatrix)
zero(::Type{DualRootSpaceElem}, ::RootSystem)
```

```@docs
root_system(::DualRootSpaceElem)
```

Basic arithmetic operations like `zero`, `+`, `-`, `*` (with rational scalars), and `==` are supported.

```@docs
coeff(::DualRootSpaceElem, ::Int)
coefficients(::DualRootSpaceElem)
```

```@docs
height(::DualRootSpaceElem)
iszero(::DualRootSpaceElem)
```

### Coroot testing
```@docs
is_coroot(::DualRootSpaceElem)
is_coroot_with_index(::DualRootSpaceElem)
is_simple_coroot(::DualRootSpaceElem)
is_simple_coroot_with_index(::DualRootSpaceElem)
is_positive_coroot(::DualRootSpaceElem)
is_positive_coroot_with_index(::DualRootSpaceElem)
is_negative_coroot(::DualRootSpaceElem)
is_negative_coroot_with_index(::DualRootSpaceElem)
```
