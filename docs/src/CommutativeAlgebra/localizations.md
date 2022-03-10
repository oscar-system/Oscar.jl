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
    Localization(S::AbsMultSet)
```
This method must be implemented with a dispatch depending on 
the concrete type of `S`.

For any concrete instance of type `AbsLocalizedRing`
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

**Note:** Any concrete type for localized rings is also required to implement 
the general [Ring Interface](@ref) of Oscar! This has not been done to a full extent 
for the previous two examples, but for `MPolyLocalizedRing`; see below.

 
### Homomorphisms for localized rings

Homomorphisms from localized rings to arbitrary algebras are of type 
```@docs
    AbsLocalizedRingHom
```
Note that, in order to be well-defined, we must have 
that ``\phi'(u) \in S`` must be a unit for every element ``u \in U``. 

The getters associated to this type which need to be implemented are 
```@docs
    domain(f::AbsLocalizedRingHom) 
    codomain(f::AbsLocalizedRingHom) 
    restricted_map(f::AbsLocalizedRingHom) 
```
Any concrete instance `f` of `AbsLocalizedRingHom` can then be applied to elements 
`a` of `domain(f)` by calling `f(a)`. 

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
   ideal(W::AbsLocalizedRing{RingType, RingElemType, MultSetType}, f::AbsLocalizedRingElem{RingType, RingElemType, MultSetType}) where {RingType, RingElemType, MultSetType}
   ideal(W::AbsLocalizedRing{RingType, RingElemType, MultSetType}, f::RingElemType) where {RingType, RingElemType, MultSetType}
   ideal(W::AbsLocalizedRing{RingType, RingElemType, MultSetType}, v::Vector{AbsLocalizedRingElem{RingType, RingElemType, MultSetType}}) where {RingType, RingElemType, MultSetType}
   ideal(W::AbsLocalizedRing{RingType, RingElemType, MultSetType}, v::Vector{RingElemType}) where {RingType, RingElemType, MultSetType}
```
for a single and a list of generators from both the `base_ring` of `W` and from `W` itself.

The minimal functionality which should be implemented for ideals is the test 
for ideal membership
```
Base.in(
    f::AbsLocalizedRingElem{RingType, RingElemType, MultSetType}, 
    I::AbsLocalizedIdeal{RingType, RingElemType, MultSetType}
  ) where {RingType, RingElemType, MultSetType}
```
and again the same for elements `f` of type `RingElemType`.

Basic operations on ideals which are already implemented on the generic level are 
```
Base.:*(I::T, J::T) where {T<:AbsLocalizedIdeal}
Base.:+(I::T, J::T) where {T<:AbsLocalizedIdeal}
```
Everything else, such as e.g. intersections of ideals, has to be implemented for the specific 
types by the user.

## Localizations of multivariate polynomial rings

Various primitive types of multiplicative sets are available, such as 
```@docs
MPolyComplementOfPrimeIdeal{
    BaseRingType, 
    BaseRingElemType,
    RingType,
    RingElemType
  } 
MPolyComplementOfKPointIdeal{
    BaseRingType,
    BaseRingElemType, 
    RingType,
    RingElemType
  } 
MPolyPowersOfElement{
    BaseRingType,
    BaseRingElemType, 
    RingType,
    RingElemType
  } 
```
Moreover, such types can be combined to products: 

**Definition (Products of multiplicative sets):**
Let ``T`` and ``U`` be multiplicative sets in a commutative ring ``R``. The product 
of ``T`` and ``U`` is defined as 
```math
  T\cdot U = \left\{ f\cdot g : f \in T \textnormal{ and }g \in U \right\}.
```
A product of multiplicative sets ``U = U_1 \cdot \dots \cdot U_r`` is called *interreduced* 
if neither one of the factors ``U_i`` is contained in one of the others ``U_j, j \neq i``.

Note that any product of multiplicative sets may be replaced by 
an interreduced one. However, such an interreduced multiplicative set is 
not unique as the following example shows:

**Example (interreduction of products of multiplicative sets):**
An interreduced factorization of a product of multiplicative sets may 
not be unique: Consider the ring ``\mathbb Z[x]`` and the multiplicative sets 
```math
  T  = \left\{(5x)^k : k \in \mathbb N_0\right\}, \quad
  T' = \left\{ x^k : k \in \mathbb N_0\right\},\quad
  S  = \left\{ c_0 \cdot x^0 : c_0 \notin 7 \mathbb Z\right\}.
```
Then ``T\cdot S = \left\{ a⋅x^k : a \notin 7\mathbb Z, k \in \mathbb N_0 \right\} = T'\cdot S``.

**Upshot:** Whenever a product is taken, some interreduced form of the 
entire product is returned. Besides the obvious simplification in 
case all factors are contained in a single one, it is difficult to 
determine which interreduction is the best one. 
Localizations of multivariate polynomial rings are of type 

The type for storing general products is 
```@docs
MPolyProductOfMultSets{
  BaseRingType,
  BaseRingElemType, 
  RingType,
  RingElemType
}
```
and products can be taken using the usual arithmetic
```@docs
product(T::AbsMPolyMultSet, U::AbsMPolyMultSet)
```
**Note:** The methods of this function naturally attempt 
to return a primitive type of multiplicative sets whenever possible. 
Hence, they are not type-stable. 


Localizations of polynomial rings are of type
```@docs
MPolyLocalizedRing{
    BaseRingType,
    BaseRingElemType,
    RingType,
    RingElemType,
    MultSetType
  } 
```
with elements of type
```@docs
MPolyLocalizedRingElem{
    BaseRingType, 
    BaseRingElemType,
    RingType,
    RingElemType, 
    MultSetType
  }
```

Ideals in localized polynomial rings are of type 
```@docs
MPolyLocalizedIdeal{BRT, BRET, RT, RET, MST}
```
Recall (see e.g. [Eis95]) that 
if ``\mathbb k`` is a Noetherian ring, any localization ``W = R[U^{-1}]`` of a 
multivariate polynomial ring ``R = \mathbb k[x_1,\dots,x_n]`` is again Noetherian and 
any ideal ``I \subset W`` is of the form ``I = I'\cdot W`` for some ideal ``I' \subset R``. 
This correspondence is not 1:1 but for any ideal ``I \subset W`` we always 
have that 
```math
  J = \left\{ x\in R : \exists u \in U : u\cdot x \in I \right\}
```
is the unique element which is maximal among all ideals ``I'`` in ``R`` for 
which ``I = I'\cdot W``. We call this the *saturated ideal* of the localization 
and it can be obtained using 
```@docs
saturated_ideal(I::MPolyLocalizedIdeal)
```
Groebner bases for the saturated ideal can be used to bring the numerators 
of any fraction ``\frac{a}{b} \in R[S^{-1}]`` into normal form and check for 
ideal membership and/or equality of elements modulo ideals in ``R[S^{-1}]``.
But for some cases, e.g. when using local orderings for localizations at 
``\mathbb k``-points, it is desirable, to have the groebner- and standard 
basis functionality available directly in the localized ring.  
To this end we have 
```@docs
LocalizedBiPolyArray{BRT, BRET, RT, RET, MST}
```
which has a monomial ordering and a `Singular` ring associated to it. 

**Note:** Transfering an element ``\frac{a}{b} \in R[S^{-1}]`` of a localized 
ring to the `Singular`-side drops all denominators and only 
the numerators appear as polynomials in `Singular`! 
Hence, a `LocalizedBiPolyArray` is not really a 1:1-correspondence 
of elements and in particular, the `Oscar` fractions can not be recovered 
from the `Singular` side. 

This is also the type returned by any Groebner- or standard basis 
computation. We make the following convention: 

**Definition:** Let ``\mathbb k[x_1,\dots,x_n][S^{-1}]`` be 
a localized polynomial ring with ``R = \mathbb k[x_1,\dots,x_n]``. 
A monomial ordering ``\geq`` is *compatible* with the localization, if 
every unit in the localization ``R_{\geq}`` is also a unit in ``R[S^{-1}]``. 

For an ideal ``I \subset R[S^{-1}]`` and a compatible monomial ordering ``\geq`` 
we say that a set of elements ``\frac{g_1}{1},\dots,\frac{g_r}{1} \in R[S^{-1}]`` is a 
groebner/standard basis for ``I`` if the elements ``g_1,\dots,g_r`` are 
a standard basis for the saturated ideal ``J`` of ``I`` in ``R``. 


**Note:** When localizing at ``\mathbb k``-points ``a = (a_1,\dots,a_n) \in \mathbb k^n``
outside the origin, the transfer of polynomials from the `Oscar` to the 
`Singular` side in `LocalizedBiPolyArray` shifts the coordinates such that 
``a`` becomes zero. A monomial ordering is always considered after application 
of such shifts. 


Groebner and standard bases of ideals can be computed for explicit 
orderings using 
```
    groebner_basis(I::MPolyLocalizedIdeal, ord::Symbol)
```
Note that depending on the type parameters of `I`, this method 
is dispatched differently, which will also lead to different 
interpretations of the ordering. For instance, for multiplicative 
sets of type `MPolyComplementOfKPointIdeal`, a shift of variables 
taking the geometric point to the origin is applied to all polynomials 
when passing to the singular side. 

If the second argument is omitted, a default ordering 
will be chosen, depending on the type of the multiplicative set. 

**Remark:** Why bother introducing Groebner and standard basis for 
localized ideals in the first place and not only work with the 
saturated ideal? The main reason is that for localizations at 
``\mathbb k``-points, the computation of the saturated ideal is 
quite expensive: It involves a primary decomposition using a 
global ordering and discarding components outside the point 
at which has been localized. Using local orderings, on the other hand, 
we can decide ideal membership or equality of elements without 
computing the saturated ideal explicitly.

The following method might be of practical interest: 
```@docs
as_affine_algebra(
  L::MPolyLocalizedRing{BRT, BRET, RT, RET, 
  MPolyPowersOfElement{BRT, BRET, RT, RET}}; 
  inverse_name::String="θ"
) where {BRT, BRET, RT, RET}
```

## Localizations of affine algebras

Let ``R = 𝕜[x₁,…,xₘ]`` be a polynomial ring, ``I ⊂ R`` some ideal 
and ``P = R/I`` its quotient. Then ``P`` is naturally an ``R``-module 
and localization of ``P`` as a ring coincides with localization 
as an ``R``-module in the sense that for every multiplicative 
set ``T ⊂ R`` there is a commutative diagram 
```math
\begin{matrix}
        R   & → & P = R/I\\
        ↓ & &       ↓ \\
  W = R[T⁻¹] & → & P[T⁻¹].
\end{matrix}
```
Observe that, moreover, for every multiplicative set 
``T' ⊂ P`` the preimage ``T`` of ``T'`` in ``R`` is also a multiplicative set. 

We may therefore treat localizations of polynomial algebras 
as localizations of modules over free polynomial rings:
and apply the following 

**Convention:** For localizations of affine algebras 
``L = (𝕜[x₁,…,xₙ]/I)[S⁻¹]`` 

  * ideals in ``L`` are given by ideals in ``W = 𝕜[x₁,…,xₙ][S⁻¹]`` containing ``I\cdot S^{-1}``.
  * the available multiplicative sets for ``L`` are exclusively those for ``𝕜[x₁,…,xₙ]``.

Note that this leads to the following differences compared to the 
standard usage of the localization interface:

 * The `base_ring` returns neither ``P``, nor ``W``, but ``R``.
 * The `BaseRingType` is the type of ``R`` and similar for 
   the other ring-based type parameters.

This is to make the data structure most accessible for 
the computational backends.

 * The type returned by `numerator` and `denominator` 
   on an element of type `MPolyQuoLocalizedRingElem` is 
   not `RingElemType`, but the type of ``P``. 

This is to comply with the purely mathematical viewpoint
where elements of localized rings are fractions of 
residue classes rather than residue classes of fractions. 


Localizations of affine algebras are realized by 
```@docs
MPolyQuoLocalizedRing{
  BaseRingType,
  BaseRingElemType,
  RingType,
  RingElemType,
  MultSetType <: AbsMultSet{RingType, RingElemType}
}
```
which have the additional methods 
```@docs
quotient_ring(L::MPolyQuoLocalizedRing)
localized_ring(L::MPolyQuoLocalizedRing)
```
returning the other two computationally important rings in its construction.

Elements of such rings are of the form 
```@docs
MPolyQuoLocalizedRingElem{
  BaseRingType, 
  BaseRingElemType,
  RingType,
  RingElemType, 
  MultSetType
}
```
In contrast to ordinary elements of a localized ring, they 
have the additional methods 
```@docs
lifted_numerator(a::MPolyQuoLocalizedRingElem)
lifted_denominator(a::MPolyQuoLocalizedRingElem)
fraction(a::MPolyQuoLocalizedRingElem)
```

