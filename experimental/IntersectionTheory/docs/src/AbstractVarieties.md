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
abstract_flag_variety(dims::Int...; base::Ring = QQ, symbol::String = "c")
```

### New Varieties From Given Varieties/Bundles

```@docs
complete_intersection(X::AbstractVariety, degs::Int...)
```

```@docs
abstract_projective_bundle(F::AbstractBundle; symbol::String = "z")
```

```@docs
abstract_hirzebruch_surface(n::Int)
```

```@docs
abstract_flag_variety(F::AbstractBundle, dims::Int...; symbol::String = "c")
```

```@docs
zero_locus_section(F::AbstractBundle; class::Bool = false)
```

```@docs
degeneracy_locus(F::AbstractBundle, G::AbstractBundle, k::Int; class::Bool=false)
```

!!! note
    Products and blowups are described elsewhere.

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
tangent_bundle(X::AbstractVariety)
```

```@docs
tautological_bundles(X::AbstractVariety)
```

```@docs
structure_map(X::AbstractVariety)
```

## Further Data Associated to an Abstract Variety


```@docs
trivial_line_bundle(X::AbstractVariety)
```

```@docs
line_bundle(X::AbstractVariety, n::RingElement)
```

```@docs
cotangent_bundle(X::AbstractVariety)
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

```@docs
basis(X::AbstractVariety)
```

```@docs
intersection_matrix(X::AbstractVariety)
```

```@docs
dual_basis(X::AbstractVariety)
```

!!! note
    If `X` is of type `AbstractVariety`, entering `total_chern_class(X)` returns the total Chern class of the tangent bundle of `X`. Similarly for entering `euler(X)`, `chern_class(X, k)`,  `todd_class(X)`, `total_pontryagin_class(X)`, `pontryagin_class(X, k)`

## Operations on Abstract Varieties

```@docs
product(X::AbstractVariety, Y::AbstractVariety)
```

!!! note
    Blowups are described in their own section.

## Integrate Chow Ring Elements

```@julia
integral(x::MPolyDecRingElem)
```

Given an element `x` of the Chow ring of an abstract variety `X`, say, return the integral of `x`.

!!! note
    If `X` has a (unique) point class, the integral will be a number (that is, a `QQFieldElem` or a function field element). Otherwise, the highests degree part of $x$ is returned (geometrically, this is the 0-dimensional part of $x$).

###### Examples

```jldoctest
julia> G = abstract_grassmannian(2, 4)
AbstractVariety of dim 4

julia> Q = tautological_bundles(G)[2]
AbstractBundle of rank 2 on AbstractVariety of dim 4

julia> E = symmetric_power(Q, 3)
AbstractBundle of rank 4 on AbstractVariety of dim 4

julia> integral(top_chern_class(E))
27

```
