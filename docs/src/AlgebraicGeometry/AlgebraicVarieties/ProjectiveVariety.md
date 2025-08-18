```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Projective Varieties
A projective variety over an algebraically closed field
is an irreducible projective algebraic set. See [Projective Algebraic Sets](@ref).

In practice we work over non-closed fields. To be called a variety
an algebraic set $V$ must stay irreducible when viewed over the algebraic closure.

In OSCAR projective varieties are [Projective schemes](@ref) and more formally defined as follows.

```@docs
AbsProjectiveVariety
```

## Constructors
```@docs
variety(I::MPolyIdeal{<:MPolyDecRingElem}; check::Bool=true)
variety(X::AbsProjectiveScheme{<:Field}; check::Bool=true)
variety(R::Ring; check::Bool=true)
variety(f::MPolyDecRingElem; check=true)
```

## Attributes
In addition to what is inherited from [Projective Algebraic Sets](@ref) and [Projective schemes](@ref), we currently have:

```@docs
sectional_genus(X::AbsProjectiveVariety)
```

## Properties
In addition to what is inherited from [Projective Algebraic Sets](@ref) and [Projective schemes](@ref), we currently have:

```@docs
is_linearly_normal(X::AbsProjectiveVariety)
```

## Methods
In addition to what is inherited from [Projective Algebraic Sets](@ref) and [Projective schemes](@ref), we currently have:

```@docs
canonical_bundle(X::AbsProjectiveVariety)
```
