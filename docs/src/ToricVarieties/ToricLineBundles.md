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

Toric line bundles `l1` and `l2` (on the same toric variety) can be tensored
by `l1*l2`. We support the `k`-th tensor power by `l1^k`. `k` can be either
an integer or valued in fmpz. The inverse of `l1` is computed by `inv(l1)`.


### Equality

Equality of toric line bundles `l1` and `l2` (on the same toric variety) is
implemented by `l1 == l2`.


## Properties

To check if a line bundle `l` is trivial, one can invoke `istrivial(l)`. Beyond this,
we support the following properties of toric line bundles:
```@docs
is_basepoint_free(l::ToricLineBundle)
isample(l::ToricLineBundle)
is_very_ample(l::ToricLineBundle)
```


## Attributes

```@docs
degree(l::ToricLineBundle)
divisor_class(l::ToricLineBundle)
toric_divisor(l::ToricLineBundle)
toric_variety(l::ToricLineBundle)
```


### Special line bundles

```@docs
anticanonical_bundle(v::AbstractNormalToricVariety)
canonical_bundle(v::AbstractNormalToricVariety)
structure_sheaf(v::AbstractNormalToricVariety)
```


## Methods

```@docs
basis_of_global_sections_via_rational_functions(l::ToricLineBundle)
basis_of_global_sections_via_homogeneous_component(l::ToricLineBundle)
```
