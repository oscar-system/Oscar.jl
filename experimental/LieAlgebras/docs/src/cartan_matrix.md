```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Cartan Matrices

Cartan matrices can be constructed from a Cartan type, and are represented as a square `ZZMatrix`.

Many functions taking a Cartan matrix as input (like [`root_system`](@ref) and [`weyl_group`](@ref)) will also accept a Cartan type as input. Both cases are semantically equivalent, but the latter may be more efficient.

!!! note
    The convention for Cartan matrices in OSCAR is $(a_{ij}) = (\langle \alpha_i^\vee, \alpha_j \rangle)$ for simple roots $\alpha_i$.

## Table of contents

```@contents
Pages = ["cartan_matrix.md"]
Depth = 2:5
```

## Constructors

```@docs
cartan_matrix(::Symbol, ::Int)
cartan_matrix(::Vector{Tuple{Symbol,Int}})
```


## Properties

```@docs
is_cartan_matrix(::ZZMatrix)
cartan_symmetrizer(::ZZMatrix)
cartan_bilinear_form(::ZZMatrix)
```


## Cartan types

The following function is used to verify the validity of a Cartan type, and is thus used in input sanitization.
```@docs
is_cartan_type(::Symbol, ::Int)
```

Given a Cartan matrix, the following functions can be used to determine its Cartan type.
```@docs
cartan_type(::ZZMatrix)
cartan_type_with_ordering(::ZZMatrix)
```
