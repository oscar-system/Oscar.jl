```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```
# Overview

The standard interfaces for groups are implemented, including `has_gens`, `number_of_generators`, `gens`, `gen`, `is_finite` and `order`.
Basic operations for group elements are also implemented. See [Basics of Groups](@ref basics_of_groups) for more details.

!!! warning
    Currently only simply connected linear algebraic groups of type A_n are supported. In the future, the functionality could be extended to support arbitrary root data.

!!! warning
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

## Special subgroups
### Root subgroups
```@docs
root_subgroup(LAG::LinearAlgebraicGroup, alpha::RootSpaceElem)
```

### Maximal torus
```@docs
maximal_torus(LAG::LinearAlgebraicGroup)
torus_element(LAG::LinearAlgebraicGroup, diag::Vector)
apply_root_to_torus_element(alpha::RootSpaceElem, t::LinearAlgebraicGroupElem)
```

### Borel subgroup
```@docs
borel_subgroup(LAG::LinearAlgebraicGroup)
```

## Bruhat decomposition
```@docs
bruhat_cell_representative(LAG::LinearAlgebraicGroup, w::WeylGroupElem)
bruhat_cell(LAG::LinearAlgebraicGroup, w::WeylGroupElem)
bruhat_decomposition(LAG::LinearAlgebraicGroup)
```

## Technicalities

