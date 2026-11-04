# Invariants of Tori

In this section, with notation as in the [introduction](@ref invariant_theory) to this chapter, $T = (K^{\ast})^m$ will be a torus of rank $m$ over an algebraically closed field $K$, acting diagonally on $V$.
As for [linearly reductive groups](@ref "How Linearly Reductive Groups and Their Representations are Given"), $T$ and its representation are defined over a subfield $k$ of $K$ supported by OSCAR, and all computations are performed over $k$.

In principle, invariants of diagonal torus actions can be computed with Derksen's algorithm from the [previous section](@ref "Invariants of Linearly Reductive Groups").
OSCAR instead uses Algorithm 4.3.1 in [DK15](@cite), which relies on polyhedral geometry and, as the authors put it, "is simpler and more efficient".

## Creating Invariant Rings

### How Tori and Their Representations are Given

```@docs
torus_group(k::Field, m::Int)
```

```@docs
rank(T::TorusGroup)
```

```@docs
base_ring(T::TorusGroup)
```

```@docs
representation_from_weights(T::TorusGroup, W::Union{ZZMatrix, Matrix{<:IntegerUnion}, Vector{<:IntegerUnion}})
```

```@docs
weights(r::RepresentationTorusGroup)
```

```@docs
group(r::RepresentationTorusGroup)
```

If `r` is a representation of a torus on $V$, then `degree(r)` and `vector_space_dim(r)` refer
to the dimension $n$ of $V$.

### Constructor for Invariant Rings

```@docs
invariant_ring(r::RepresentationTorusGroup)
```

## Basic Data Associated to Invariant Rings

If `RT` is the invariant ring $k[x_1, \dots, x_n]^T$ of a torus $T$, then

- `group(RT)` refers to $T$,
- `coefficient_ring(RT)` to $k$, and
- `polynomial_ring(RT)` to $k[x_1, \dots, x_n]$.

## Fundamental Systems of Invariants

```@docs
fundamental_invariants(RT::TorGroupInvarRing)
```

## Invariant Rings as Affine Algebras

```@docs
affine_algebra(RT::TorGroupInvarRing)
```
