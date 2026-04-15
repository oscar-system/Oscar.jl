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

### Pullback and pushforward semantics

The **pullback** $f^*$ is a ring homomorphism: it preserves products and sends $1$ to $1$. It is
always available once a map is constructed.

The **pushforward** $f_*$ is only a group homomorphism (it does *not* preserve products in general),
but it satisfies the **projection formula**:

$$f_*(f^*(y) \cdot x) = y \cdot f_*(x) \quad \text{for all } x \in \mathrm{N}^*(X),\; y \in \mathrm{N}^*(Y).$$

When the pushforward is not specified explicitly, OSCAR attempts to compute it from the pullback
using the projection formula together with the intersection pairing on $Y$. This succeeds
automatically when:

- **Y** is a point or a curve;
- all algebraic classes in $\mathrm{N}^*(Y)_{\mathbb Q}$ are known (e.g., projective spaces, Grassmannians, complete intersections);
- the flag `:alg` is set on $Y$.

In other cases, a warning is issued and the result may be incorrect.

```jldoctest
julia> P2 = abstract_projective_space(2);

julia> P5 = abstract_projective_space(5, symbol = "H");

julia> h = gens(P2)[1]
h

julia> i = map(P2, P5, [2*h]); # Veronese embedding

julia> pushforward(i, h) # pushforward works via projection formula
2*H^4

julia> pushforward(i, P2(1)) # pushforward of the fundamental class
4*H^3

```

### Inclusions and `extend_inclusion`

For an inclusion $i\colon X \hookrightarrow Y$, where the class $[X]$ is not already representable
in $\mathrm{N}^*(Y)_{\mathbb Q}$, the pushforward $i_*$ cannot be expressed in the existing Chow
ring of $Y$. Setting `inclusion = true` in the `map` constructor calls `extend_inclusion` to create a
modified variety $Y^+$ with extra generators added to the Chow ring so that $i_*$ becomes well-defined.
The structure map $Y^+ \to Y$ gives access to the enlarged ring while preserving the original
ring via pullback.

## Types

The OSCAR type for abstract variety maps is `AbstractVarietyMap`.

## Constructors

```@docs
map(X::AbstractVariety, Y::AbstractVariety, f_pullback::Vector, f_pushforward = nothing; inclusion::Bool = false, symbol::String = "x")
```

```@docs
extend_inclusion(i::AbstractVarietyMap; symbol::String = "e")
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

