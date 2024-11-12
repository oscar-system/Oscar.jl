```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Cartan Matrices

Cartan matrices can be constructed from a Cartan type, and are represented as a square `ZZMatrix`.

Many functions taking a Cartan matrix as input (like [`root_system`](@ref) and [`weyl_group`](@ref)) will also accept a Cartan type as input. Both cases are semantically equivalent, but the latter may be more efficient.

!!! note
    The convention for Cartan matrices in OSCAR is $(a_{ij}) = (\langle \alpha_i^\vee, \alpha_j \rangle)$ for simple roots $\alpha_i$.

```@docs
cartan_matrix(::Symbol, ::Int)
cartan_matrix(::Vector{Tuple{Symbol,Int}})
```

```@docs
is_cartan_matrix(::ZZMatrix)
cartan_symmetrizer(::ZZMatrix)
cartan_bilinear_form(::ZZMatrix)
```

```@docs
cartan_type(::ZZMatrix)
cartan_type_with_ordering(::ZZMatrix)
```

```@docs
is_cartan_type(::Symbol, ::Int)
```
