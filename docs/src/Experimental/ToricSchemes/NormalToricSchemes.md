```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["NormalToricSchemes.md"]
```

# Normal Toric Schemes

## Constructors

We can construct a toric scheme as follows:

```julia
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> toric_scheme = ToricCoveredScheme(P2)
Scheme of a toric variety with fan spanned by RayVector{fmpq}[[1, 0], [0, 1], [-1, -1]]
```


## Properties

We currently support the following properties:

```@docs
is_smooth(X::ToricCoveredScheme)
```


## Attributes

We currently support the following attributes:

```@docs
underlying_scheme(X::ToricCoveredScheme)
normal_toric_variety(X::ToricCoveredScheme)
fan(X::ToricCoveredScheme)
```
