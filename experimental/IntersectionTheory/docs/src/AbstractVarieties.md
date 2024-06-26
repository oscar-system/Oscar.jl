```@meta
CurrentModule = Oscar
```

# Abstract Varieties

## Types

AbsVariety <: Variety

## Constructors

```@docs
abstract_variety(n::Int, A::MPolyDecRingOrQuo)
```

```@docs
abstract_point()
```

### Specialized Constructors

```@docs
abstract_projective_space(n::Int; base::Ring = QQ, symbol::String = "h")
```

```@docs
abstract_grassmannian(k::Int, n::Int; bott::Bool = false, weights = :int, base::Ring = QQ, symbol::String = "c")
```

```@docs 
abstract_flag_variety(dims::Int...; bott::Bool = false)
```

For examples of constructors starting from already given abstract varieties, vector bundles, or homomorphisms
between vector bundles see  [`product`](@ref), [`complete_intersection`](@ref), [`abstract_projective_bundle`](@ref), 
[`section_zero_locus`](@ref), and [`degeneracy_locus`](@ref).

## Underlying Data of an Abstract Variety

An abstract variety is made up from (a selection of) the data discussed here:

```@docs
dim(X::AbstractVariety)
```

```@docs
chow_ring(X::AbstractVariety)
```

```@docs
base(X::AbstractVariety)
```

```@docs
point_class(X::AbstractVariety)
```

```@docs
trivial_line_bundle(X::AbstractVariety)
```

```@docs
tangent_bundle(X::AbstractVariety)
```

```@docs
tautological_bundles(X::AbstractVariety)
```

```@docs
structure_map(X::AbstractVariety)
```

## Further Data Associated to Abstract Varieties

```@docs
cotangent_bundle(X::AbstractVariety)
```

```@docs
line_bundle(X::AbstractVariety, n::RingElement)
```

```@docs
canonical_class(X::AbstractVariety)
```

```@docs
canonical_bundle(X::AbstractVariety)
```

```@docs
degree(X::AbstractVariety)
```

If `X` is of type `AbstractVariety` or `TnVariety`, entering `total_chern_class(X)` returns the total Chern class of the tangent bundle of `X`.
Similarly for entering `euler(X)`, `chern_class(X, k)`,  `todd_class(X)`, `total_pontryagin_class(X)`, `pontryagin_class(X, k)`

## Operations on Abstract Varieties

```@docs
product(X::AbstractVariety, Y::AbstractVariety)
```
