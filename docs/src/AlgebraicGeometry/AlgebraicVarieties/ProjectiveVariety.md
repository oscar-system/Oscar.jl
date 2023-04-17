```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

# Projective Varieties
A projective variety over an algebraically closed field
is an irreducible projective algebraic set. See [Projective Algebraic Sets](@ref).

In practice we work over non-closed fields. To be called a variety
an algebraic set $V$ must stay irreducible when viewed over the algebraic closure.

In Oscar projective varieties are [Projective schemes](@ref) and more formally defined as follows.

```@docs
AbsProjectiveVariety
```

## Constructors
```@docs
projective_variety(I::MPolyIdeal{<:MPolyDecRingElem}; check::Bool=true)
projective_variety(X::AbsProjectiveScheme{<:Field}; check::Bool=true)
projective_variety(R::Ring; check::Bool=true)
projective_variety(f::MPolyDecRingElem; check=true)
```

## Attributes
So far all are inherited from [Projective Algebraic Sets](@ref) and [Projective schemes](@ref).

## Properties
So far all are inherited from [Projective Algebraic Sets](@ref) and [Projective schemes](@ref).

## Methods
So far all are inherited from [Projective Algebraic Sets](@ref) and [Projective schemes](@ref).
