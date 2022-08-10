```@meta
CurrentModule = Oscar
```

```@contents
Pages = ["ToricLineBundles.md"]
```


# Toric Line Bundles


## Constructors

### Generic constructors

```@docs
ToricLineBundle(v::AbstractNormalToricVariety, c::Vector{T}) where {T <: IntegerUnion}
ToricLineBundle(v::AbstractNormalToricVariety, d::ToricDivisor)
```

### Tensor products

Toric line bundles can be tensored via `*`. The `n`-th tensor power can be computed via `^n`.
In particular, `^(-1)` computes the inverse of a line bundle. Alternatively, one can compute
the inverse by invoking `inv`.

### Special line bundles

```@docs
anticanonical_bundle(v::AbstractNormalToricVariety)
canonical_bundle(v::AbstractNormalToricVariety)
structure_sheaf(v::AbstractNormalToricVariety)
```


## Properties

Equality of toric line bundles can be tested via `==`.

To check if a toric line bundle is trivial, one can invoke `is_trivial`. Beyond this,
we support the following properties of toric line bundles:
```@docs
is_basepoint_free(l::ToricLineBundle)
is_ample(l::ToricLineBundle)
is_very_ample(l::ToricLineBundle)
```


## Attributes

```@docs
degree(l::ToricLineBundle)
divisor_class(l::ToricLineBundle)
toric_divisor(l::ToricLineBundle)
toric_variety(l::ToricLineBundle)
```


## Methods

```@docs
basis_of_global_sections_via_rational_functions(l::ToricLineBundle)
basis_of_global_sections_via_homogeneous_component(l::ToricLineBundle)
```
