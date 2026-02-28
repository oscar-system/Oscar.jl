```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```
# Overview

The standard interfaces for groups are implemented, including `has_gens`, `number_of_generators`, `gens`, `gen`, `is_finite` and `order`.
Basic operations for group elements are also implemented. See [Basics of Groups](@ref basics_of_groups) for more details.

!!! warning
    Currently only root systems of type A_n are supported. In the future the functionality could be extended to support arbitrary root data.

We are computing everything using the standard example (i.e. ``SL_n``) as a matrix group with standard torus, standard Borel etc.
Most functionality is currently limited to finite fields.

## Table of contents
```@contents
Pages = ["overview.md"]
Depth = 2:5
```

## Constructing groups
```@docs
linear_algebraic_group(type::Symbol, n::Int, k::Field)
linear_algebraic_group(rs::RootSystem, k::Field)
```

## Constructing elements
```@docs
linear_algebraic_group_elem(LAG::LinearAlgebraicGroup{C}, MGE::MatGroupElem{C}) where {C<:FieldElem}
linear_algebraic_group_elem(LAG::LinearAlgebraicGroup{C}, m::MatElem{C}) where {C<:FieldElem}
```

## Root subgroups
```@docs
root_subgroup_generator(LAG::LinearAlgebraicGroup, alpha::RootSpaceElem)
```

```@docs
root_subgroup(LAG::LinearAlgebraicGroup, alpha::RootSpaceElem)
```

## Tori
We work with the standard torus of diagonal elements.
```@docs
maximal_torus(LAG::LinearAlgebraicGroup)
```

```@docs
torus_element(LAG::LinearAlgebraicGroup, diag::Vector)
```

```@docs
apply_root_to_torus_element(alpha::RootSpaceElem, t::LinearAlgebraicGroupElem)
```

## Bruhat decomposition
We work with the standard Borel of upper triangular matrices.
```@docs
representative_of_root_in_group(LAG::LinearAlgebraicGroup, alpha::RootSpaceElem)
```

```@docs
borel_subgroup(LAG::LinearAlgebraicGroup)
```

```@docs
bruhat_cell_representative(LAG::LinearAlgebraicGroup, w::WeylGroupElem)
```

```@docs
bruhat_cell(LAG::LinearAlgebraicGroup, w::WeylGroupElem)
```

```@docs
bruhat_decomposition(LAG::LinearAlgebraicGroup)
```

