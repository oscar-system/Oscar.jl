```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```


# Morphisms of affine schemes


## Constructors

### General constructors

```@docs
morphism(X::AbsAffineScheme, Y::AbsAffineScheme, f::Vector{<:RingElem}; check::Bool=true)
```

### Special constructors 

```@docs
id_hom(X::AbsAffineScheme{<:Any, <:MPolyRing})
inclusion_morphism(X::AbsAffineScheme, Y::AbsAffineScheme; check::Bool=true)
compose(f::AbsAffineSchemeMor, g::AbsAffineSchemeMor)
restrict(f::AffineSchemeMor, U::AbsAffineScheme, V::AbsAffineScheme)
```


## Attributes

### General attributes

```@docs
domain(f::AbsAffineSchemeMor)
codomain(f::AbsAffineSchemeMor)
pullback(f::AbsAffineSchemeMor)
graph(f::AbsAffineSchemeMor)
```

### Special attributes

In addition to the standard getters and methods for instances
of `AffineSchemeMor`, we also have
```@docs
image_ideal(f::ClosedEmbedding)
```

### Undocumented

The following functions do exist but are currently undocumented:
- `underlying_morphism`,
- `complement_ideal`,
- `complement_scheme`,
- `preimage`,
- `inverse`,
- various type getters.


## Properties

```@docs
is_isomorphism(f::AbsAffineSchemeMor)
is_inverse_of(f::AbsAffineSchemeMor, g::AbsAffineSchemeMor)
is_identity_map(f::AbsAffineSchemeMor)
```


## Methods

```@docs
fiber_product(f::AbsAffineSchemeMor, g::AbsAffineSchemeMor)
product(X::AbsAffineScheme, Y::AbsAffineScheme)
simplify(X::AbsAffineScheme{<:AbstractAlgebra.Field})
```
