```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

# Affine Varieties
An affine variety over an algebraically closed field
is an irreducible affine algebraic set. See [Affine Algebraic Sets](@ref).

In practice we work over non-closed fields. To be called a variety
an algebraic set $V$ must stay irreducible when viewed over the algebraic closure.

In Oscar varieties are [Affine schemes](@ref) and more formally defined as follows.
```@docs
AbsAffineVariety
```

## Constructors
```@docs
affine_variety(I::MPolyIdeal; check=true)
affine_variety(X::Spec{<:Field}; check::Bool=true)
affine_variety(R::MPolyAnyRing; check=true)
```

## Attributes
So far all are inherited from [Affine Algebraic Sets](@ref) and [Affine schemes](@ref).

## Properties
So far all are inherited from [Affine Algebraic Sets](@ref) and [Affine schemes](@ref).

## Methods
So far all are inherited from [Affine Algebraic Sets](@ref) and [Affine schemes](@ref).

