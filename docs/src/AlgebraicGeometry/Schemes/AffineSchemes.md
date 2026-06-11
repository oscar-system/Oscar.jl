```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```


# Affine schemes

Let ``\mathbb k`` be a commutative noetherian base ring
(in practice: an algebraic extension of ``\mathbb Q`` or ``\mathbb F_p``).
We support functionality for affine schemes ``X = \mathrm{Spec}(R)`` over ``\mathbb k``.
Currently, we support rings ``R`` of type `MPolyRing`, `MPolyQuoRing`,
`MPolyLocRing`, and `MPolyQuoLocRing`
defined over the integers, a finite field or algebraic field extensions of ``\mathbb Q``


## Constructors

### General constructors

Besides `spec(R)` for `R` of either one of the types `MPolyRing`, `MPolyQuoRing`, `MPolyLocRing`, or
`MPolyQuoLocRing`, we have the following constructors:
```@docs
spec(R::MPolyRing, I::MPolyIdeal)
spec(R::MPolyRing, U::AbsMPolyMultSet)
spec(R::MPolyRing, I::MPolyIdeal, U::AbsMPolyMultSet)
```
See [`inclusion_morphism(::AbsAffineScheme, ::AbsAffineScheme)`](@ref) for a way to obtain the ideal ``I`` from ``X = \mathrm{Spec}(R, I)``.

### Affine n-space

```@docs
affine_space(kk::BRT, n::Int; variable_name=:x) where {BRT<:Ring}
affine_space(kk::BRT, var_names::AbstractVector{<:VarName}) where {BRT<:Ring}
```

### Closed subschemes

```@docs
subscheme(X::AbsAffineScheme, f::Vector{<:RingElem})
subscheme(X::AbsAffineScheme, I::Ideal)
```

### Intersections

```@docs
Base.intersect(X::AbsAffineScheme{BRT, <:Ring}, Y::AbsAffineScheme{BRT, <:Ring}) where {BRT<:Ring}
```

### Open subschemes

```@docs
hypersurface_complement(X::AbsAffineScheme, f::RingElem)
hypersurface_complement(X::AbsAffineScheme, f::Vector{<:RingElem})
```

### Closure

```@docs
closure(X::AbsAffineScheme, Y::AbsAffineScheme)
```


## Attributes

### Ambient affine space

Most affine schemes in OSCAR ``X = \mathrm{Spec}(R)``
over a ring ``B``, come with an embedding into an
affine space ``\mathbb{A}_B``.
More precisely, `ambient_space(X)` is defined for `X = spec(R)` if `R`
is constructed from a polynomial ring.
In particular ``\mathrm{Spec}(\mathbb{Z})`` or ``\mathrm{Spec}(\mathbb{k})`` for ``\mathbb k``
a field do not have an ambient affine space.

```@docs
ambient_space(X::AbsAffineScheme)
```

### Other attributes

```@docs
base_ring(X::AbsAffineScheme)
codim(X::AbsAffineScheme)
ambient_embedding(X::AbsAffineScheme)
dim(X::AbsAffineScheme)
name(X::AbsAffineScheme)
OO(X::AbsAffineScheme)
```

### Type getters

We support functions which return the types of
schemes, associated rings, and their elements. See the
source code for details.


## Properties

```@docs
is_open_embedding(X::AbsAffineScheme, Y::AbsAffineScheme)
is_closed_embedding(X::AbsAffineScheme, Y::AbsAffineScheme)
isempty(X::AbsAffineScheme)
is_subscheme(X::AbsAffineScheme, Y::AbsAffineScheme)
```


## Methods
```@docs
tangent_space(X::AbsAffineScheme{<:Field}, P::AbsAffineRationalPoint)
is_normal(X::AbsAffineScheme; check::Bool=true)
normalization(X::AbsAffineScheme; check::Bool=true, algorithm=:equidimDec)
```

### Comparison

Two schemes ``X`` and ``Y`` can be compared if their ambient affine spaces are equal.
In particular ``X`` and ``Y`` are considered equal (`==`)
if and only if the identity morphism of their ambient affine space induces an
isomorphism of ``X`` and ``Y``.
For ``X`` and ``Y`` with different ambient affine space `X==Y` is always `false`.

### Auxiliary methods

```@docs
is_non_zero_divisor(f::RingElem, X::AbsAffineScheme)
```
