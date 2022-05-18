```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

# Localizations of commutative rings

Suppose ``R`` is a commutative ring with unit and ``U \subset R`` is a
*multiplicatively closed subset,* that is,
```math
s, t \in U \;\Rightarrow \; s\cdot t \in U \;\text{ and }\; 1 \in U.
```
Then we can form the *localization* of ``R`` at ``U,``
```math
    R[U^{-1}] = \left\{ \frac{p}{q} : p,q \in R, \, q \in U \right\},
```
which comes equipped with its standard arithmetic for fractions. See, for instance, [Eis95](@cite).

Oscar provides a general  localization framework, originally intended to be used 
with multivariate polynomial rings ``R`` over some base field ``\mathbb k``, but also 
applicable to more general commutative rings.

In the case of  a polynomial ring ``R``, the localization framework provides
the structure to integrate algorithms which, using Gröbner bases, allow one to solve  ideal- and module-theoretic
questions concerning ``R[U^{-1}]`` by computations over ``R``. This makes it important to regard localizations
``R[U^{-1}]`` as rings with a history of creation from the original pair ``U \subset R``. 

!!! note
    In the local context, following Hironaka and Grauert, it is also common to use the name
    *standard basis* instead of Gröbner basis. See, for example, [GP08](@cite).

## The Localization Interface

For the convenience of developers, we describe the interface that needs to be implemented in order
to create a concrete instance of localized rings. At the same time, we illustrate the use of the interface
by examples featuring localizations of multivariate polynomial rings, a concrete instance already
realized in OSCAR. Functionality for this instance, which includes functions such as
`MPolyComplementOfPrimeIdeal` for creating multiplicatively closed subsets, will be
discussed in a subsequent section. 

### Localized Rings

Multiplicatively closed subsets are derived from the abstract type
```@docs
    AbsMultSet{RingType, RingElemType}
```
The basic functionality that has to be implemented for any concrete type derived from 
this is to be able to check containment of elements in multiplicatively closed subsets via
```@docs
    in(f::RingElemType, U::AbsMultSet{RingType, RingElemType}) where {RingType, RingElemType}
```
This is supposed to be an extension of the methods of the function `Base.in`.

Any concrete type of localized rings should then be a subtype of
```@docs
    AbsLocalizedRing{RingType, RingElemType, MultSetType}
```
The basic way to construct localized rings is to first 
specify a multiplicative set `U`, and then call 
```@docs
    Localization(U::AbsMultSet)
```
This method must be implemented with a dispatch depending on 
the concrete type of `U`.

For any concrete instance of type `AbsLocalizedRing`,
the following methods must be implemented:
```@docs
    base_ring(Rloc::AbsLocalizedRing) 
    inverted_set(Rloc::AbsLocalizedRing)
```

### Elements of Localized Rings

Also, conversion of fractions to elements of localized rings must be implemented in the form 
`(W::AbsLocalizedRing{RingType, RingElemType, MultSetType})(a::RingElemType) where {RingType, RingElemType, MultSetType}`, taking ``a`` to the element ``\frac{a}{1}``.
For more general fractions one needs
`(W::AbsLocalizedRing{RingType, RingElemType, MultSetType})(a::RingElemType, b::RingElemType) where {RingType, RingElemType, MultSetType}`, mapping a pair ``(a, b)`` to the fraction ``\frac{a}{b}``.


The *elements* of localized rings must be derived from 
```@docs
    AbsLocalizedRingElem{RingType, RingElemType, MultSetType}
```
For any concrete instance `f` of `AbsLocalizedRingElem` there must be the following 
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

 
### Homomorphisms From Localized Rings

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

### Ideals in Localized Rings

Finitely generated ideals in localizations should be of type
```@docs
    AbsLocalizedIdeal{LocRingElemType<:AbsLocalizedRingElem}
```
The required getter methods are
```julia
    gens(I::AbsLocalizedIdeal)
    base_ring(I::AbsLocalizedIdeal)
```
The constructors to be implemented are
```
   ideal(W::AbsLocalizedRing, f::AbsLocalizedRingElem) 
   ideal(W::AbsLocalizedRing, v::Vector{LocalizedRingElemType}) where {LocalizedRingElemType<:AbsLocalizedRingElem}
```

The minimal functionality which must be implemented for ideals is the following
for ideal membership and writing elements as linear combinations of the generators:
```
    Base.in(f::AbsLocalizedRingElem, I::AbsLocalizedIdeal)
    coordinates(f::AbsLocalizedRingElem, I::AbsLocalizedIdeal)
```

## Localizations of Multivariate Rings

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
```julia
    MPolyLocalizedIdeal{
        LocRingType<:MPolyLocalizedRing, 
        LocRingElemType<:MPolyLocalizedRingElem
      } <: AbsLocalizedIdeal{LocRingElemType}
```
Recall (see e.g. [Eis95]) that 
if ``\mathbb k`` is a Noetherian ring, any localization ``W = R[U^{-1}]`` of a 
multivariate polynomial ring ``R = \mathbb k[x_1,\dots,x_n]`` is again Noetherian and 
any ideal ``I \subset W`` is of the form ``I = J\cdot W`` for some ideal ``J \subset R``. 
There is an ambiguity in the choice of such ``J``, but we always have that
for any choice of ``J`` as above
```math
  J:U = \left\{ x\in R : \exists u \in U : u\cdot x \in J \right\} = \left\{ x \in R : x//1 \in I \right\}
```
is the unique element which is maximal among all ideals ``J'`` in ``R`` for 
which ``I = J'\cdot W``. We call this the *saturated ideal* ``J' = J:U`` of ``I``
and it can be obtained using 
```@docs
saturated_ideal(I::MPolyLocalizedIdeal)
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

