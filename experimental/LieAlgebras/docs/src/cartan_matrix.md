```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Cartan Matrices

```@docs
cartan_matrix(::Symbol, ::Int)
cartan_matrix(::Vector{Tuple{Symbol,Int}})
cartan_matrix(::Tuple{Symbol,Int}...)
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
