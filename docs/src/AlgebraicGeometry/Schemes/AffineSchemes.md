```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["AffineSchemes.md"]
```


# Affine schemes

Let ``\mathbb k`` be a commutative noetherian base ring
(in practice: an algebraic extension of ``\mathbb Q`` or ``\mathbb F_p``).
We support functionality for affine schemes ``X = \mathrm{Spec}(R)`` over ``\mathbb k``.
Currently, we support rings ``R`` of type `MPolyRing`, `MPolyQuo`,
`MPolyLocalizedRing`, and `MPolyQuoLocalizedRing`
defined over the integers or algebraic field extensions of ``\mathbb Q``


## Constructors

### General constructors

Besides `Spec(R)` for `R` of either one of the types `MPolyRing`, `MPolyQuo`, `MPolyLocalizedRing`, or 
`MPolyQuoLocalizedRing`, we have the following constructors:
```@docs
Spec(R::MPolyRing, I::MPolyIdeal)
Spec(R::MPolyRing, U::AbsMPolyMultSet)
Spec(R::MPolyRing, I::MPolyIdeal, U::AbsMPolyMultSet)
```

### Affine n-space

```@docs
affine_space(kk::BRT, n::Int; variable_name="x") where {BRT<:Ring}
affine_space(kk::BRT, var_symbols::Vector{Symbol}) where {BRT<:Ring}
```

### Closed subschemes

```@docs
subscheme(X::AbsSpec, f::Vector{<:RingElem})
subscheme(X::AbsSpec, I::Ideal)
```

### Intersections

```@docs
Base.intersect(X::AbsSpec{BRT, <:Ring}, Y::AbsSpec{BRT, <:Ring}) where {BRT<:Ring}
```

### Open subschemes

```@docs
hypersurface_complement(X::AbsSpec, f::RingElem)
hypersurface_complement(X::AbsSpec, f::Vector{<:RingElem})
```

### Closure

```@docs
closure(X::AbsSpec, Y::AbsSpec)
```


## Attributes

### Ambient affine space

Most affine schemes in Oscar ``X = \mathrm{Spec}(R)``
over ``\mathbb k``, come with an embedding into an
affine space.
```@docs
ambient_affine_space(X::AbsSpec)
```
Some such as ``Spec(ZZ)`` or ``Spec(k)`` for ``\mathbb k``
a field do not have an ambient affine space.

### Other attributes

```@docs
base_ring(X::AbsSpec)
codim(X::AbsSpec)
ambient_closure_ideal(X::AbsSpec{<:Any, <:MPolyRing})
dim(X::AbsSpec)
name(X::Spec)
OO(X::AbsSpec)
```

### Type getters

We support functions which return the types of
schemes, associated rings, and their elements. See the 
source code for details.


## Properties

```@docs
is_open_embedding(X::AbsSpec, Y::AbsSpec)
is_closed_embedding(X::AbsSpec, Y::AbsSpec)
isempty(X::AbsSpec)
issubset(X::AbsSpec, Y::AbsSpec)
```


## Methods

### Comparison

Two schemes $X$ and $Y$ can be compared if their ambient affine spaces are equal.
In particular $X$ and $Y$ are considered equal (`==`)
if and only if the identity morphism of their ambient affine space induces an
isomorphism of ``X`` and ``Y``.
For ``X`` and ``Y`` with different ambient affine space `X==Y` is always `false`.

### Auxilliary methods

```@docs
is_non_zero_divisor(f::RingElem, X::AbsSpec)
```
