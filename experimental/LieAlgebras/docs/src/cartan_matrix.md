```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Cartan Matrices

```@docs
cartan_matrix(::Symbol, ::Int)
cartan_matrix(::Tuple{Symbol,Int}...)
is_cartan_matrix(::ZZMatrix; generalized::Bool)
cartan_symmetrizer(::ZZMatrix; check::Bool)
cartan_bilinear_form(::ZZMatrix; check::Bool)
cartan_type(::ZZMatrix; check::Bool)
cartan_type_with_ordering(::ZZMatrix; check::Bool)
```
