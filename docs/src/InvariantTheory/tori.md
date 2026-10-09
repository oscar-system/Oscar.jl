# Invariants of Tori
In this section, with notation as in the [introduction](@ref invariant_theory) to this chapter, $T =(K^{\ast})^m$ will be a torus of rank $m$
over a field $K$. To compute invariants of diagonal torus actions, OSCAR makes use of Algorithm 4.3.1  in [DK15](@cite) which,
in particular, relies on algorithmic means from polyhedral geometry.

## Creating Invariant Rings

### How Tori  and Their Representations are Given

```@docs
 torus_group(F::Field, n::Int)
```

```@docs
rank(T::TorusGroup)
```

```@docs
field(T::TorusGroup)
```

```@docs
representation_from_weights(T::TorusGroup, W::Union{ZZMatrix, Matrix{<:Integer}, Vector{<:Int}})
```

```@docs
group(r::RepresentationTorusGroup)
```

### Constructor for Invariant Rings

```@docs
invariant_ring(r::RepresentationTorusGroup)
```


## Fundamental Systems of Invariants

```@docs
fundamental_invariants(RT::TorGroupInvarRing)
```


## Invariant Rings as Affine Algebras

```@docs
affine_algebra(RT::TorGroupInvarRing)
```
