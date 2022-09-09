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

## Affine schemes over noetherian base rings

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
subscheme(X::Spec, f::RingElem)
hypersurface_complement(X::Spec, f::MPolyElem)
```
Containments can be checked using 
```@docs 
    issubset(X::Spec, Y::Spec)
    is_open_embedding(X::Spec, Y::Spec)
    is_closed_embedding(X::Spec, Y::Spec)
```
The closure of a subscheme can be computed via
```@docs 
    closure(X::Spec, Y::Spec)
```
Among the basic functionality for affine schemes we have 
```@docs
    product(X::Spec, Y::Spec)
```
