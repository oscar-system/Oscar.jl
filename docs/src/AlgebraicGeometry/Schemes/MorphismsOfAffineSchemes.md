```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["MorphismsOfAffineSchemes.md"]
```


# Morphisms of affine schemes


## Constructors

### General constructors

```@docs
SpecMor(X::AbsSpec, Y::AbsSpec, f::Vector{<:RingElem}; check::Bool=true)
```

### Special constructors (@id spec_mor_special_constructors)

```@docs
identity_map(X::AbsSpec{<:Any, <:MPolyRing})
inclusion_morphism(X::AbsSpec, Y::AbsSpec; check::Bool=true)
compose(f::AbsSpecMor, g::AbsSpecMor)
restrict(f::SpecMor, U::AbsSpec, V::AbsSpec)
```


## Attributes

### General attributes

```@docs
domain(f::AbsSpecMor)
codomain(f::AbsSpecMor)
pullback(f::AbsSpecMor)
graph(f::AbsSpecMor)
```

### Special attributes

In addition to the standard getters and methods for instances
of `SpecMor`, we also have
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
is_isomorphism(f::AbsSpecMor)
is_inverse_of(f::AbsSpecMor, g::AbsSpecMor)
is_identity_map(f::AbsSpecMor)
```


## Methods

```@docs
fiber_product(f::SpecMor{SpecType, SpecType, <:Any}, g::SpecMor{SpecType, SpecType, <:Any}) where {SpecType<:StdSpec}
product(X::AbsSpec, Y::AbsSpec)
simplify(X::AbsSpec{<:AbstractAlgebra.Field})
```
