```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

# Localizations of commutative rings

Suppose ``R`` is a commutative ring with unit and ``S \subset R`` is a *multiplicatively 
closed set* containing ``1 \in R``. Then we can form the *localization* of ``R`` at ``S``
```math
    R[S^{-1}] = \left\{ \frac{p}{q} : p,q \in R, \, q \in S \right\},
```
with its standard arithmetic for fractions. See, for instance, [Eis95] for an account on localizations.

Oscar provides a general framework for such localizations, originally intended to be used 
with multivariate polynomial rings ``R`` over some base field ``\mathbb k``, but also 
applicable to more general commutative rings.

In the case of polynomials, the localization framework provides the structure for 
certain algorithms using standard bases. Note that, in general, localizations of 
polynomial algebras are not finitely generated 
as algebras over ``\mathbb k``; for instance when localizing at some maximal 
ideal ``\mathfrak m \subset R``. However, many ideal- and module-theoretic questions in the localization 
``R[S^{-1}]``, such as e.g. the ideal membership, can be transformed to questions on 
ideals and modules over the base ring ``R`` and then solved using Groebner- or standard-basis 
techniques. This makes it important to regard localizations ``R[S^{-1}]`` as rings with 
a history of creation from the original pair ``S \subset R``. 

## The localization interface

### Localized rings

The interface that needs to be implemented for any concrete 
instance of localized rings is the following. 
Multiplicatively closed sets are derived from the abstract type
```@docs
    AbsMultSet{RingType, RingElemType}
```
The basic functionality that has to be implemented for any concrete type derived from 
this is to be able to check containment of elements via
```@docs
    in(f::RingElemType, S::AbsMultSet{RingType, RingElemType}) where {RingType, RingElemType}
```
This is supposed to be an extension of the methods of the function `Base.in`.

A localized ring should then be derived from 
```@docs
    AbsLocalizedRing{RingType, RingElemType, MultSetType}
```
The basic way to construct localized rings is to first 
specify a multiplicative set `S` and then call 
```@docs
    localize_at(S::AbsMultSet)
```
This method must be implemented with a dispatch depending on 
the concrete type of `S`.

For any concrete instance of type ```AbsLocalizedRing``` 
the following methods must be implemented:
```@docs
    base_ring(W::AbsLocalizedRing) 
    inverted_set(W::AbsLocalizedRing)
```
Also, conversion of fractions to elements of localized rings must be implemented in the form 
`(W::AbsLocalizedRing{RingType, RingElemType, MultSetType})(a::RingElemType) where {RingType, RingElemType, MultSetType}`, taking ``a`` to the element ``\frac{a}{1}``.
For more general fractions one needs
`(W::AbsLocalizedRing{RingType, RingElemType, MultSetType})(a::RingElemType, b::RingElemType) where {RingType, RingElemType, MultSetType}`, mapping a pair ``(a, b)`` to the fraction ``\frac{a}{b}``.


The *elements* of localized rings must be derived from 
```@docs
    AbsLocalizedRingElem{RingType, RingElemType, MultSetType}
```
For any concrete instance `F` of `AbsLocalizedRingElem` there must be the following 
methods:
```@docs
    numerator(f::AbsLocalizedRingElem) 
    denominator(f::AbsLocalizedRingElem) 
    parent(f::AbsLocalizedRingElem)
```
A default version of the arithmetic is implemented on the generic level using the above 
functionality and the arithmetic for the original ring. 
Depending on the actual concrete instance, one might wish to provide more fine-tuned methods, 
starting e.g. by implementing 
```@docs
    reduce_fraction(f::AbsLocalizedRingElem)
```
Note that this is called after *every* arithmetic operation (addition, multiplication,...), 
so the computations carried out here should be computationally cheap.

Two toy examples for implementations of this interface for localizations 
of the integers ``\mathbb Z`` and rings of the form ``\mathbb Z/n\mathbb Z`` 
can be found in the test files 
`test/Rings/integer-localizations.jl` and `test/Rings/nmod-localizations.jl`.
 
### Ideals in localized rings

One of the main reasons to implement localizations in the first place 
is that this process preserves the property of a ring to be Noetherian; 
which is crucial for computer algebra. In this regard, we have 
```@docs
    AbsLocalizedIdeal{RingType, RingElemType, MultSetType} 
```
The required getter methods are
```@docs
    gens(I::AbsLocalizedIdeal)
    base_ring(I::AbsLocalizedIdeal)
```
The constructors to be implemented are
```
   ideal(W::AbsLocalizedRing{RingType, RingElemType, MultSetType}, f::AbsLocalizedRingElem{RingType, RingElemType, MultSetType})
   ideal(W::AbsLocalizedRing{RingType, RingElemType, MultSetType}, f::RingElemType)
   ideal(W::AbsLocalizedRing{RingType, RingElemType, MultSetType}, v::Vector{AbsLocalizedRingElem{RingType, RingElemType, MultSetType}})
   ideal(W::AbsLocalizedRing{RingType, RingElemType, MultSetType}, v::Vector{RingElemType})
```
for a single and a list of generators from both the `base_ring` of `W` and from `W` itself.

The minimal functionality which should be implemented for ideals is the test 
for ideal membership
```@docs
Base.in(
    f::AbsLocalizedRingElem{RingType, RingElemType, MultSetType}, 
    I::AbsLocalizedIdeal{RingType, RingElemType, MultSetType}
  ) where {RingType, RingElemType, MultSetType}
```
and again the same for elements `f` of type `RingElemType`.

Basic operations on ideals which are already implemented on the generic level are 
```@docs
Base.:*(I::T, J::T) where {T<:AbsLocalizedIdeal}
Base.:+(I::T, J::T) where {T<:AbsLocalizedIdeal}
```
Everything else, such as e.g. intersections of ideals, has to be implemented for the specific 
types by the user.
