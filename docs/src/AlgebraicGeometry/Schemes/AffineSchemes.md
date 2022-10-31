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
We support functionality for affine schemes ``X = Spec(R)`` over ``\mathbb k``.
Currently, we support rings ``R`` of type `MPolyRing`, `MPolyQuo`,
`MPolyLocalizedRing`, and `MPolyQuoLocalizedRing`
defined over the integers or algebraic field extensions of ``\mathbb Q``


## Constructors

### General constructors

```@docs
Spec(R::MPolyRing, I::MPolyIdeal)
Spec(R::MPolyRing, U::AbsMPolyMultSet)
Spec(R::MPolyRing, I::MPolyIdeal, U::AbsMPolyMultSet)
```

### Copy constructors

```@docs
Spec(X::Spec)
```

### Affine n-space

```@docs
affine_space(kk::BRT, n::Int; variable_name="x") where {BRT<:Ring}
affine_space(kk::BRT, var_symbols::Vector{Symbol}) where {BRT<:Ring}
```

### Standard spectrum

```@docs
standard_spec(X::AbsSpec{<:Any, <:MPolyRing})
standard_spec(X::AbsSpec{<:Any, <:MPolyQuo})
standard_spec(X::AbsSpec{<:Any, <:MPolyLocalizedRing})
standard_spec(X::AbsSpec{<:Any, <:MPolyQuoLocalizedRing})
```

### Closed subschemes

```@docs
subscheme(X::AbsSpec, f::Vector{<:RingElem})
subscheme(X::AbsSpec, I::Ideal)
```

### Intersections

```@docs
Base.intersect(X::AbsSpec{BRT, <:MPolyRing}, Y::AbsSpec{BRT, <:MPolyRing}) where {BRT<:Ring}
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

### Ambient ring

For most affine schemes ``X = \mathrm{Spec}(R)``
over ``\mathbb k``, there is a 'governing' polynomial
``\mathbb k``-algebra ``P = \mathbb{k}[x_1,\dots,x_n]``
in the following sense:
```@docs
    ambient_ring(X::AbsSpec)
```
For instance, this is the case whenever ``R`` is a quotient
ring of ``P``, a localization of ``P``, or a localization
of a quotient ring of ``P``; but also for
power series rings in multiple variables.

### Other attributes

```@docs
base_ring(X::AbsSpec)
codim(X::AbsSpec)
defining_ideal(X::AbsSpec{<:Any, <:MPolyRing})
dim(X::AbsSpec)
name(X::Spec)
OO(X::AbsSpec)
strict_modulus(X::AbsSpec)
```

### Type getters

We support functions which return the types of
schemes and associated rings.


## Properties

```@docs
is_open_embedding(X::AbsSpec, Y::AbsSpec)
is_closed_embedding(X::AbsSpec, Y::AbsSpec)
isempty(X::AbsSpec)
issubset(X::AbsSpec, Y::AbsSpec)
```


## Methods

### Comparison

Schemes can be compared based on their `ambient_ring`.

Two subschemes `X` and `Y` of the same ambient affine space `A` are considered equal
if and only if the identity morphism of `A` induces an isomorphism of `X` and `Y`.
This is implemented by `X==Y`. For `X` and `Y` in different ambient affine spaces,
`X==Y` is always `false`.

### Auxilliary methods

```@docs
is_non_zero_divisor(f::RingElem, X::AbsSpec{<:Ring, <:MPolyRing})
```
