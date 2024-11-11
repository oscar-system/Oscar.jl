```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Root systems

```@docs
root_system(::ZZMatrix)
root_system(::Symbol, ::Int64)
root_system(::Vector{Tuple{Symbol, Int64}})
```

```@docs
cartan_matrix(::RootSystem)
weyl_group(::RootSystem)
```

```@docs
is_simple(::RootSystem)
rank(::RootSystem)
weyl_vector(::RootSystem)
```

```@docs
has_root_system_type(::RootSystem)
root_system_type(::RootSystem)
root_system_type_with_ordering(::RootSystem)
```

```@docs
fundamental_weight(::RootSystem, ::Int64)
fundamental_weights(::RootSystem)
```

```@docs
number_of_roots(::RootSystem)
number_of_positive_roots(::RootSystem)
number_of_simple_roots(::RootSystem)
```

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


## Root space elements

```@docs
RootSpaceElem(::RootSystem, ::Vector{<:RationalUnion})
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

```@docs
reflect(::RootSpaceElem, ::Int)
reflect!(::RootSpaceElem, ::Int)
```


## Dual root space elements

```@docs
DualRootSpaceElem(::RootSystem, ::Vector{<:RationalUnion})
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

```@docs
is_coroot(::RootSpaceElem)
is_coroot_with_index(::RootSpaceElem)
is_simple_coroot(::RootSpaceElem)
is_simple_coroot_with_index(::RootSpaceElem)
is_positive_coroot(::RootSpaceElem)
is_positive_coroot_with_index(::RootSpaceElem)
is_negative_coroot(::RootSpaceElem)
is_negative_coroot_with_index(::RootSpaceElem)
```

## Weight lattice elements

```@docs
WeightLatticeElem(::RootSystem, ::Vector{<:IntegerUnion})
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
```

```@docs
reflect(::WeightLatticeElem, ::Int)
reflect!(::WeightLatticeElem, ::Int)
```

```@docs
conjugate_dominant_weight(::WeightLatticeElem)
conjugate_dominant_weight_with_left_elem(::WeightLatticeElem)
conjugate_dominant_weight_with_right_elem(::WeightLatticeElem)
```
