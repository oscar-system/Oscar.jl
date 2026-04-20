```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Abstract bundles

An *abstract bundle* on an abstract variety $X$ is determined its Chern character (or, equivalently, by its rank and total Chern class).
Abstract bundles support the standard operations on vector bundles in algebraic geometry: direct sum, tensor product, duals,
determinant bundles, exterior and symmetric powers, as well as pullback and pushforward along abstract variety maps.
They also carry the usual characteristic classes: Chern classes, Segre classes, Todd class, and Pontryagin classes.

The arithmetic operations `+`, `-`, `*` on abstract bundles correspond to direct sum, formal difference,
and tensor product, respectively. In particular, multiplying a bundle by an integer `n` gives the
direct sum of `n` copies.

### Grothendieck ring and virtual bundles

Abstract bundles live in the Grothendieck ring $\mathrm{K}^0(X)$ of vector bundles on $X$.
In this ring, every element can be written as a formal difference $[E] - [F]$ of genuine
bundles. The Chern character $\mathrm{ch}\colon\mathrm{K}^0(X) \to \mathrm{N}^*(X)_{\mathbb Q}$ is
a ring homomorphism that maps the Grothendieck ring (after tensoring with $\mathbb Q$)
onto the Chow ring.

In OSCAR, an `AbstractBundle` is determined by its Chern character, so virtual bundles
with zero (and negative) rank are fully supported:

```jldoctest
julia> P2 = abstract_projective_space(2);

julia> T = tangent_bundle(P2);

julia> F = T - 2*OO(P2); # a virtual bundle of rank 0

julia> rank(F)
0

julia> chern_character(F)
3//2*h^2 + 3*h

```

### Segre classes

The *total Segre class* $s(E) = c(E)^{-1}$ is the formal inverse of the total Chern class
and arises naturally as the fundamental class of a projective bundle. For a bundle $E$ of
rank $r$ the individual Segre classes satisfy $s_k(E) = (-1)^k c_k(E^\vee)$ when $k \le r$,
but the relation $s(E)\cdot c(E) = 1$ also determines the higher Segre classes uniquely.

```jldoctest
julia> P3 = abstract_projective_space(3);

julia> T = tangent_bundle(P3);

julia> total_segre_class(T)
-20*h^3 + 10*h^2 - 4*h + 1

julia> total_segre_class(T) * total_chern_class(T) # must be 1
1

```

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
