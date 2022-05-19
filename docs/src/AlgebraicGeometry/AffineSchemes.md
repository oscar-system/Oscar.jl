```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

# General schemes

Arbitrary schemes over a base ring ``\mathbb k`` which are given by means 
of their affine patches and glueings are instances of 
```@docs
Scheme{BaseRingType<:Ring, BaseRingElemType<:RingElement}
```

# Affine schemes over noetherian base rings

Let ``\mathbb k`` be a commutative noetherian base ring 
(in practice: an algebraic extension of ``\mathbb Q`` or ``\mathbb F_p``). 
An affine scheme ``X`` over ``\mathbb k`` is given as 
```math
    X = \mathrm{Spec} \left(\mathbb k[x_1,\dots,x_n]/I\right)
```
for some polynomial ring ``R = \mathbb k[x_1,\dots,x_n]`` and an ideal ``I \subset R``.
More generally, we can form the `Spec` of any localized affine algebra ``L`` as an 
instance of 
```@docs
Spec{BRT, BRET, RT, RET, MST}
```
Internally, this merely stores an instance ``L`` of `MPolyQuoLocalizedRing`.
This ring can be obtained using 
```@docs
OO(X::Spec)
```
One of the main reasons to allow such general schemes is that *principal open subsets* 
``U = \{ f \neq 0\}`` of affine schemes ``X``, ``f \in \mathcal O_X(X)`` are again 
affine; basically due to Rabinowitsch's trick:
```math
    U = \mathrm{Spec} \left(\mathcal O_X(X)[f^{-1}]\right) = \mathrm{Spec} \left(\mathcal O_X(X)[t]/\langle 1 - t\cdot f \rangle\right).
```
This flexibility allows us to introduce and handle subschemes using the natural 
relations on their rings. Thus, we can, for instance, define subschemes via 
```@docs
subscheme(X::Spec{BRT, BRET, RT, RET, MST}, f::RET) where {BRT, BRET, RT, RET, MST}
hypersurface_complement(X::Spec{BRT, BRET, RT, RET, MST}, f::RET) where {BRT, BRET, RT, RET, MST<:MPolyPowersOfElement{BRT, BRET, RT, RET}}
```
Containments can be checked using 
```@docs 
issubset(
  X::Spec{BRT, BRET, RT, RET, MST1}, 
  Y::Spec{BRT, BRET, RT, RET, MST2}
) where {BRT, BRET, RT, RET, MST1<:MPolyPowersOfElement{BRT, BRET, RT, RET}, MST2<:MPolyPowersOfElement{BRT, BRET, RT, RET}}
is_closed_embedding(
  X::Spec{BRT, BRET, RT, RET, MST1}, 
  Y::Spec{BRT, BRET, RT, RET, MST2}
) where {BRT, BRET, RT, RET, MST1<:MPolyPowersOfElement{BRT, BRET, RT, RET}, MST2<:MPolyPowersOfElement{BRT, BRET, RT, RET}}
is_open_embedding(
  X::Spec{BRT, BRET, RT, RET, MST1}, 
  Y::Spec{BRT, BRET, RT, RET, MST2}
) where {BRT, BRET, RT, RET, MST1<:MPolyPowersOfElement{BRT, BRET, RT, RET}, MST2<:MPolyPowersOfElement{BRT, BRET, RT, RET}}
```
The closure of a subscheme can be computed via
```@docs 
closure(
  X::Spec{BRT, BRET, RT, RET, MST1}, 
  Y::Spec{BRT, BRET, RT, RET, MST2}
) where {BRT, BRET, RT, RET, MST1<:MPolyPowersOfElement{BRT, BRET, RT, RET}, MST2<:MPolyPowersOfElement{BRT, BRET, RT, RET}}
```
Among the basic functionality for affine schemes we have 
```@docs
product(X::Spec{BRT, BRET, RT, RET, MST}, Y::Spec{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST<:MPolyPowersOfElement}
```

# Projective Schemes

In Oscar we provide the functionality to handle *relative projective space* and their subschemes 
over various base rings as presented in [Har77](@cite).
```@docs
    ProjectiveScheme{CoeffRingType, CoeffRingElemType, RingType, RingElemType}
```
Currently, one can take the ring ``A`` to be either a computable field such as `QQ` or `GF(p)`, 
a multivariate polynomial ring over the latter, or, with a view towards the affine schemes 
above, an instance of `MPolyQuoLocalizedRing`, the type of rings used for the structure 
sheaves. This ring of coefficients can be asked for using 
```@docs
    base_ring(X::ProjectiveScheme)
```
and the various other getters are 
```@docs
    fiber_dimension(X::ProjectiveScheme)
    homogeneous_coordinate_ring(X::ProjectiveScheme)
    generators_of_defining_ideal(X::ProjectiveScheme)
```
Internally, projective varieties are modeled via their affine cones. 
```@docs
    affine_cone(X::ProjectiveScheme) 
```
In order to be able to fluently switch between the homogeneous coordinate 
ring ``S`` of a projective scheme ``X`` and the coordinate ring ``R`` of 
its affine cone ``C(X)``, we provide the following methods:
```@docs
    homogeneous_coordinates(X::ProjectiveScheme)
    convert_to_fraction(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::RET) where {CRT<:MPolyQuoLocalizedRing, CRET, RT, RET}
    convert_to_homog_polys(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::MPolyLocalizedRingElem) where {CRT<:MPolyRing, CRET, RT, RET}
```

  
