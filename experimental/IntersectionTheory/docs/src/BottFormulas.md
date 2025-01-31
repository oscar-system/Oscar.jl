```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Bott's Formula

## Abstract Varieties With a Torus Action

### Types

`TnVariety`

### Constructors

```@docs
tn_variety(n::Int, points::Vector{Pair{P, Int}}) where P
```

### Specialized Constructors

```@docs
tn_grassmannian(k::Int, n::Int; weights = :int)
```

```@docs
tn_flag_variety(dims::Int...; weights = :int)
```

### Underlying Data of an Abstract Variety With a Torus Action

`dimension(X::TnVariety)`
`points(X::TnVariety)`
`tangent_bundle(X::TnVariety)`
`bundles(X::TnVariety)`

### Further Data Associated to an Abstract Variety With a Torus Action

!!! note
    If `X` is of type `TnVariety`, entering `total_chern_class(X)` returns the total Chern class of the tangent bundle of `X`.
    Similarly for entering `chern_class(X, k)`,
	
`trivial_line_bundle(X::TnVariety)`
`euler_number(X::TnVariety)`
`cotangent_bundle(X::TnVariety)`
`euler_number(X::TnVariety)`

`total_chern_class(X::TnVariety)`
`chern(k::Int, X::TnVariety)`

## Equivariant Abstract Bundles Under a Torus Action

### Types

### Constructors

### Specialized Constructors

### Underlying Data of an Equivariant Abstract Bundle

### Further Data Associated to an Equivariant Abstract Bundle

## Integrating ...



## Example: Linear Subspaces on Hypersurfaces

```@docs
linear_subspaces_on_hypersurface(k::Int, d::Int; bott::Bool = true)
```
