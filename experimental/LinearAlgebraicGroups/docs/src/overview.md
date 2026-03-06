```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```
# Overview

The standard interfaces for groups are implemented, including `has_gens`, `number_of_generators`, `gens`, `gen`, `is_finite` and `order`.
Basic operations for group elements are also implemented. See [Basics of Groups](@ref basics_of_groups) for more details.
Additionally we provide `normalizer`, `centralizer`, `is_conjugate`, `conjugate_group` and `in`.

!!! warning
    Currently only simply connected linear algebraic groups of type ``A_n`` are supported. In the future, the functionality could be extended to support arbitrary root data.

!!! warning
    Most functionality is currently limited to finite fields.

## Table of contents
```@contents
Pages = ["overview.md"]
Depth = 2:5
```

## Constructing groups
```@docs
linear_algebraic_group(rs::RootSystem, k::Field)
```

## Constructing elements
```@docs
linear_algebraic_group_elem(LAG::LinearAlgebraicGroup{C}, MGE::MatGroupElem{C}) where {C<:FieldElem}
```

## Properties
Use `root_system` to obtain the root system of the linear algebraic group.
`base_ring` returns the field the group is defined over.

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
The used types are `LinearAlgebraicGroup{C}` for the groups and `LinearAlgebraicGroupElem{C}` for their elements.
The type parameter `C` is the element type of the base field.
These types wrap a `MatGroup{C}` or a `MatGroupElem{C}` respectively, but can
store additional information about the group that would not be possible if we
would just use the existing matrix groups.
