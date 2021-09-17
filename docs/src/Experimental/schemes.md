```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["schemes.md"]
```

# Introduction

This package is supposed to provide the basic functionality for 
schemes in OSCAR.
This comprises affine schemes $U = R/I$ with 
$R = k[x_1,\dots,x_n]$ a polynomial algebra over a base ring $k$ (usually a field), 
as well as 
schemes $X = \bigcup_{i=1}^N U_i$ which arise from glueings of 
affine patches $U_i$. In addition, we provide the basic objects such as 
line and vector bundles, coherent sheaves, and the standard operations 
on them. 

The guiding idea in the development of this package is to make the handling of
projective schemes $X$ independent of the choice of an embedding $X \hookrightarrow 
\mathbb P^n$ and its description via some homogeneous coordinate ring. 
In practice, the dimension $n$ of the ambient space for such embeddings is so high 
that most Groebner basis driven algorithms applied to X will not finish within 
a human life span. By glueing $X = \bigcup_{i=1}^N U_i$ from local affine patches we aim to 
keep the number of variables $n(i)$, which are necessary for the description of the patches 
$U_i = k[x_1,\dots,x_{n(i)}]/J_i$, reasonably low so that algorithms involving Groebner 
basis computations on those patches will finish in reasonable time. 
For the computation of global invariants of a covered scheme $X$, however, we need 
new tools which may then resort to the outputs of the local computations.


# Schemes and their morphisms

## Affine schemes 

The very basic building block is 
```@docs 
Spec{S,T,U}
```
```@docs
   base_ring(::Spec)
```

```@docs
   ambient_ring(::Spec)
```

```@docs
   defining_ideal(::Spec)
```

Spec is a subtype of
```@docs
AffineScheme{S, T <: MPolyRing, U <: MPolyElem}
```
which is itself derived from 
```@docs
Scheme{ S <: Ring }
```

Another important instance of affine schemes is 
```@docs
    SpecPrincipalOpen{S,T,U}
```
```@docs
   parent(::SpecPrincipalOpen)
```
```@docs
   denom(::SpecPrincipalOpen)
```
```@docs
   denoms(::SpecPrincipalOpen{S,T,U}) where {S <: Ring, T <: MPolyRing, U<:MPolyElem}
```

The explicit description of a SpecPrincipalOpen D via an affine algebra 
can be asked for by the user via 
```@docs
    ambient_ring(::SpecPrincipalOpen)
```
and 
```@docs
    defining_ideal(::SpecPrincipalOpen)
```
To really derive an instance of Spec from this, simply call 
```@docs
    Spec(::SpecPrincipalOpen)
```
Given the rings above, one can ask for the inverses of 
the localized elements via
```@docs
   inverted_element(::SpecPrincipalOpen)
```

## Morphisms of affine schemes

```@docs
    AffSchMorphism{S,Tdom, Udom, Tcod, Ucod}
```

```@docs
   domain(::AffSchMorphism)
```

```@docs
   codomain(::AffSchMorphism)
```

```@docs
   imgs_frac(::AffSchMorphism)
```

```@docs
   pullback(::AffSchMorphism)
```

## Covered schemes

```@docs
    CoveredScheme
```

A covering itself relies on 

```@docs
   Glueing
```

These are then organized in 

```@docs 
    Covering
```

# Vector bundles and coherent sheaves

```@docs
   LineBundle
```
Such a bundle is described by an "antisymmetric" matrix of transition 
functions for the specified covering. 
```@docs
   covering(::LineBundle)
```
```@docs
   transitions(::LineBundle)
```
Note that not every line bundle on an affine 
scheme is necessarily trivial! In general one will 
need to even cover an affine scheme by finer patches 
in order to represent it as a line bundle of this form.

With basically the same implementation one has 
```@docs
    VectorBundle
```
which we separated from the line bundles because of the 
special role that the latter play for tensor multiplication.
For vector bundles, the transition functions take values 
in the space of matrices:
```@docs
    transitions(::VectorBundle)
```

## Basic operations on vector bundles 
```@docs
    direct_sum(::VectorBundle,::VectorBundle)
```
This functionality has been extended in the obvious way to also 
work for line bundles. 

# Constructing new schemes from old

```@docs
    projectivize(::VectorBundle)
```


