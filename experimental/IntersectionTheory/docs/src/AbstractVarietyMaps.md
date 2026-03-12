```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Abstract variety maps

An *abstract variety map* $f\colon X \to Y$ encodes a morphism between abstract varieties at the level
of intersection theory. It is determined by:

- the **pullback** $f^*\colon \mathrm{N}^*(Y)_{\mathbb Q}\to \mathrm{N}^*(X)_{\mathbb Q}$, a ring homomorphism on Chow rings;
- optionally, a **pushforward** $f_*\colon \mathrm{N}^*(X)_{\mathbb Q}\to \mathrm{N}^*(Y)_{\mathbb Q}$, a group homomorphism satisfying the projection formula.

When the pushforward is not given explicitly, it can sometimes be computed automatically via the
projection formula, provided enough information about the Chow ring of $Y$ is available (see the
documentation of `map` for details).

Abstract variety maps also carry a *relative tangent bundle* $\mathrm{T}_f$, which satisfies
$\mathrm{T}_X = f^* \mathrm{T}_Y \oplus \mathrm{T}_f$, and may carry a *relative polarization* $\mathcal{O}_f(1)$.

## Types

The OSCAR type for abstract variety maps is `AbstractVarietyMap`.

## Constructors

```@docs
map(X::AbstractVariety, Y::AbstractVariety, f_pullback::Vector, f_pushforward = nothing; inclusion::Bool = false, symbol::String = "x")
```

```@docs
identity_map(X::AbstractVariety)
```

## Underlying data of an abstract variety map

An abstract variety map is made up from (a selection of) the data discussed here:

```@docs
domain(f:: AbstractVarietyMap)
```

```@docs
codomain(f:: AbstractVarietyMap)
```

```@docs
dim(f::AbstractVarietyMap)
```

```@docs
pullback(f::AbstractVarietyMap, y::MPolyDecRingOrQuoElem)
```

```@docs
pushforward(f::AbstractVarietyMap, x::MPolyDecRingOrQuoElem)
```

```@docs
tangent_bundle(f::AbstractVarietyMap)
```

## Further data associated to an abstract variety map

```@docs
cotangent_bundle(f::AbstractVarietyMap)
```

```@docs
todd_class(f::AbstractVarietyMap)
```

## Operations on abstract variety maps

```@docs
compose(f::AbstractVarietyMap, g::AbstractVarietyMap)
```

