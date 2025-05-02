```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
```


# Toric Line Bundles


## Constructors

### Generic constructors

```@docs
toric_line_bundle(v::NormalToricVarietyType, picard_class::FinGenAbGroupElem)
toric_line_bundle(v::NormalToricVarietyType, picard_class::Vector{T}) where {T <: IntegerUnion}
toric_line_bundle(v::NormalToricVarietyType, d::ToricDivisor)
toric_line_bundle(d::ToricDivisor)
toric_line_bundle(v::NormalToricVarietyType, dc::ToricDivisorClass)
toric_line_bundle(dc::ToricDivisorClass)
```

### Tensor products

Toric line bundles can be tensored via `*`. The `n`-th tensor power can be computed via `^n`.
In particular, `^(-1)` computes the inverse of a line bundle. Alternatively, one can compute
the inverse by invoking `inv`.

### Special line bundles

```@docs
anticanonical_bundle(v::NormalToricVarietyType)
canonical_bundle(v::NormalToricVarietyType)
structure_sheaf(v::NormalToricVarietyType)
trivial_line_bundle(v::NormalToricVarietyType)
```


## Properties

Equality of toric line bundles can be tested via `==`.

To check if a toric line bundle is trivial, one can invoke `is_trivial`. Beyond this,
we support the following properties of toric line bundles:
```@docs
is_ample(l::ToricLineBundle)
is_basepoint_free(l::ToricLineBundle)
is_immaculate(l::ToricLineBundle)
is_very_ample(l::ToricLineBundle)
```


## Attributes

```@docs
degree(l::ToricLineBundle)
picard_class(l::ToricLineBundle)
toric_divisor(l::ToricLineBundle)
toric_divisor_class(l::ToricLineBundle)
toric_variety(l::ToricLineBundle)
```


## Methods

```@docs
basis_of_global_sections_via_rational_functions(l::ToricLineBundle)
basis_of_global_sections_via_homogeneous_component(l::ToricLineBundle)
generic_section(l::ToricLineBundle)
```
