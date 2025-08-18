```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Affine Varieties
An affine variety is an algebraic set such that $X(K)$ is irreducible for $k \subseteq K$ an algebraic closure.
See [Affine Algebraic Sets](@ref).

In OSCAR varieties are implemented as special instances of [Affine schemes](@ref) and more formally defined as follows.
```@docs
AbsAffineVariety
```
Functionality which is not (yet) provided by a variety-specific implementation, falls back to the appropriate functionality of schemes.

## Constructors
```@docs
variety(I::MPolyIdeal; check=true)
variety(X::AbsAffineScheme{<:Field}; is_reduced=false, check::Bool=true)
variety(R::MPolyAnyRing; check=true)
```

## Attributes
So far all are inherited from [Affine Algebraic Sets](@ref) and [Affine schemes](@ref).

## Properties
So far all are inherited from [Affine Algebraic Sets](@ref) and [Affine schemes](@ref).

## Methods
So far all are inherited from [Affine Algebraic Sets](@ref) and [Affine schemes](@ref).

