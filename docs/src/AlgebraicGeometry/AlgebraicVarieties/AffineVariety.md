```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

# Affine Varieties
An affine variety is a geometrically irreducible affine algebraic set $X$.
This means that $X(K)$ is irreducible for $k \subseteq K$ an algebraic closure.
See [Affine Algebraic Sets](@ref).

In Oscar varieties can are schemes [Affine schemes](@ref) and more formally defined as follows.
```@docs
AbsAffineVariety
```

## Constructors
```@docs
variety(I::MPolyIdeal; check=true)
variety(X::AbsSpec{<:Field}; is_reduced=false, check::Bool=true)
variety(R::MPolyAnyRing; check=true)
```

## Attributes
So far all are inherited from [Affine Algebraic Sets](@ref) and [Affine schemes](@ref).

## Properties
So far all are inherited from [Affine Algebraic Sets](@ref) and [Affine schemes](@ref).

## Methods
So far all are inherited from [Affine Algebraic Sets](@ref) and [Affine schemes](@ref).

