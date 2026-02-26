```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Abstract bundles

An *abstract bundle* on an abstract variety $X$ is determined by its rank and its Chern character
(or, equivalently, its total Chern class). Abstract bundles support the standard operations on vector bundles
in algebraic geometry: direct sum, tensor product, duals, determinant bundles, exterior and symmetric powers,
as well as pullback and pushforward along abstract variety maps.
They also carry the usual characteristic classes: Chern classes, Segre classes, Todd class, and Pontryagin classes.

The arithmetic operations `+`, `-`, `*` on abstract bundles correspond to direct sum, formal difference,
and tensor product, respectively. In particular, multiplying a bundle by an integer `n` gives the
direct sum of `n` copies.

## Types

The OSCAR type for abstract vector bundles is `AbstractBundle`.

## Constructors

```@docs
abstract_bundle(X::AbstractVariety, ch::Union{MPolyDecRingElem, MPolyQuoRingElem})
```

## Underlying data of an abstract bundle

An abstract bundle is made up from (a selection of) the data discussed here:

```@docs
parent(F::AbstractBundle)
```

```@docs
chern_character(F::AbstractBundle)
```

```@docs
rank(F::AbstractBundle)
```

```@docs
total_chern_class(F::AbstractBundle)
```

## Further data associated to an abstract bundle

```@docs
chern_class(F::AbstractBundle, k::Int)
```

```@docs
top_chern_class(F::AbstractBundle)
```

```@docs
total_segre_class(F::AbstractBundle)
```

```@docs
segre_class(F::AbstractBundle, k::Int)
```

```@docs
 todd_class(F::AbstractBundle)
```

```@docs
total_pontryagin_class(F::AbstractBundle)
```

```@docs
pontryagin_class(F::AbstractBundle, k::Int)
```

```@docs
euler_characteristic(F::AbstractBundle)
```

```@docs
hilbert_polynomial(F::AbstractBundle)
```

## Operations on abstract bundles

The arithmetic operations `+`, `-`, `*` on abstract bundles correspond to direct sum, formal
difference, and tensor product, respectively. Multiplying a bundle by an integer `n` gives the
direct sum of `n` copies.

```jldoctest
julia> P3 = abstract_projective_space(3);

julia> 4*OO(P3, 1) - OO(P3) == tangent_bundle(P3) # Euler sequence
true

```

```@docs
dual(F::AbstractBundle)
```

```@docs
det(F::AbstractBundle)
```

```@docs
exterior_power(F::AbstractBundle, k::Int)
```

```@docs
symmetric_power(F::AbstractBundle, k::Int)
```

```@docs
pullback(f::AbstractVarietyMap, F::AbstractBundle)
```

```@docs
pushforward(f::AbstractVarietyMap, F::AbstractBundle)
```

## Tests on abstract bundles

```@docs
==(F::AbstractBundle, G::AbstractBundle)
```
