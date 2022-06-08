```@meta
CurrentModule = Oscar 
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["localizations.md"]
```

# Localization

Let ``R`` be a commutative ring with 1, and let ``U \subset R`` be a
*multiplicatively closed subset,* that is,
```math
s, t \in U \;\Rightarrow \; s\cdot t \in U \;\text{ and }\; 1 \in U.
```
Then we can form the *localization* of ``R`` at ``U,``  that is,  the set
```math
    R[U^{-1}] = \left\{ \frac{p}{q} : p,q \in R, \, q \in U \right\},
```
equipped with its standard arithmetic for fractions. See, for instance, [Eis95](@cite).

In this section, we describe OSCAR functionality for handling localizations of multivariate polynomial
rings, as well as localizations of affine algebras, at various types of multiplicatively closed subsets. For the convenience of
the developer, in an appendix, we outline a general framework for creating new concrete instances of localized
rings in OSCAR, describing relevant abstract types and indicating the functionality to be implemented.

!!! note
    All OSCAR types discussed in this section are parametrized. For simplicity of the presentation, we omit corresponding details.


## Localizing Multivariate Rings

Most functions for handling localizations of multivariate polynomial rings rely on Gröbner bases. More precisely,
if ``R`` is a multivariate polynomial ring, Gröbner bases allow one to solve  ideal- and module-theoretic
questions concerning ``R[U^{-1}]`` by computations over ``R``.

!!! note
    Recall that in OSCAR, Gröbner bases are implemented for multivariate polynomial rings over fields (exact fields supported by OSCAR) and for multivariate polynomial rings over the integers.

!!! note
    In the local context, following Hironaka and Grauert, it is often common to use the name *standard basis* instead of Gröbner basis. See, for example, [GP08](@cite).
	

### Types

All OSCAR types for multiplicatively closed subsets of commutative rings with unit
belong to the abstract type `AbsMultSet`.  For multiplicatively
closed subsets of multivariate polynomial rings, there are the concrete descendants
`MPolyComplementOfKPointIdeal`, `MPolyComplementOfPrimeIdeal`, and `MPolyPowersOfElement`.

The general abstract type for localizations of commutative rings with unit is `AbsLocalizedRing`.
For localizations of multivariate polynomial rings, there is the concrete subtype `MPolyLocalizedRing`.

### Constructors

#### Multiplicatively Closed Subsets

In accordance with the above mentioned types, we have the following constructors
for multiplicatively closed subsets of multivariate polynomial rings.

```@docs
complement_of_ideal(R::MPolyRing, a::Vector)
```

```@docs
powers_of_element(f::MPolyElem)
```

It is also possible to build products of multiplicatively closed sets already given:

```@docs
product(T::AbsMPolyMultSet, U::AbsMPolyMultSet)
```

Containment in multiplicatively closed subsets can be checked via the `in` function as indicated below:

##### Examples

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
S = complement_of_ideal(R, [0, 0 ,0])
y in S
P = ideal(R, [x])
T = complement_of_ideal(P)
y in T
f = x
U = powers_of_element(f)
x^3 in U
(1+y)*x^2 in product(S, U)
```


#### Localized Rings

```@docs
Localization(U::AbsMPolyMultSet)
```
	
### Data associated to Localized Rings

If `Rloc` is the localization of a multivariate polynomial ring `R`  at a multiplicatively closed subset
`U` of `R`, then

- `base_ring(Rloc)` refers to `R`, and
- `inverted_set(Rloc)` to `U`.

###### Examples

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
P = ideal(R, [x])
U = MPolyComplementOfPrimeIdeal(P)
Rloc, _ = Localization(U);
R === base_ring(Rloc)
U === inverted_set(Rloc)
```

### Elements of Localized Rings

#### Types

The general abstract type for elements of localizations of commutative rings with unit is `AbsLocalizedRingElem`.
For elements of localizations of multivariate polynomial rings, there is the concrete subtype `MPolyLocalizedRingElem`.


#### Creating Elements of Localized Rings

Elements of a localized multivariate polynomial ring ``R[U^{-1}]`` are created as (fractions of) images of elements of $R$
under the localization map or by directly coercing (pairs of) elements of $R$ into fractions.

##### Examples

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
P = ideal(R, [x])
U = MPolyComplementOfPrimeIdeal(P)
Rloc, iota = Localization(U);
f = iota(x)
f == Rloc(x)
g = Rloc(y, z)
f+g
f*g
```

#### Data Associated to Elements of Localized Rings

Given an element `f` of a localized multivariate ring polynomial `Rloc`, 
- `parent(f)` refers to `Rloc`,
- `numerator(f)` to the numerator  of the internal representative of `f`, and
- `denominator(f)` to the denominator of the internal representative of `f`.

##### Examples

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
P = ideal(R, [x])
U = MPolyComplementOfPrimeIdeal(P)
Rloc, iota = Localization(U);
f = iota(x)//iota(y)
parent(f)
g = iota(y)//iota(z)
numerator(f*g)
denominator(f*g)
```

### Homomorphisms from Localized Rings

Let ``R[U^{-1}]`` be the localization of a multivariate polynomial ring ``R`` at a
multiplicatively closed subset `U` of ``R``. Then, by the universal property of localization, each
homomorphism from ``R[U^{-1}]`` to a commutative ring  $S$ with unit is determined as the extension
of a ring homomorphism $R \to S$ sending elements of $U$ to units in $S$. In OSCAR, such
homomorphisms are of type `MPolyLocalizedRingHom`. They are created using the following constructor:

```@docs
hom(Rloc::MPolyLocalizedRing, S::Ring, F::Map)
```
	
Given a ring homomorphism `PHI` from `Rloc` to `S` as above, `domain(PHI)` and `codomain(PHI)`
refer to `Rloc` and `S`, respectively. Moreover, the composition of `PHI` with the localization map
is recovered as follows:

```@docs
restricted_map(PHI::MPolyLocalizedRingHom)
```
	  
### Ideals in Localized Rings
	  
#### Types

Ideals in localized rings are of concrete type `MPolyLocalizedIdeal`.

#### Constructors

Given a localization `Rloc` of a multivariate polynomial ring `R`, and given a vector `V` of elements of
`Rloc` (of `R`),  the ideal of `Rloc` which is generated by (the images) of the entries of `V`
is created by entering `ideal(Rloc, V)`. 

#### Data Associated to Ideals

If `I` is an ideal of a localized multivariate polynomial ring  `Rloc`, then

- `base_ring(I)` refers to `Rloc`,
- `gens(I)` to the generators of `I`,
- `ngens(I)` to the number of these generators, and
- `gen(I, k)` as well as `I[k]` to the `k`-th such generator.

#### Operations on Ideals

If `I`,  `J` are ideals of a localized multivariate polynomial ring  `Rloc`, then

- `I^k` refers to the `k`-th  power of `I`, 
- `I+J`, `I*J`,  and `intersect(I, J)` to the sum, product, and intersection of `I` and  `J`, and
- `quotient(I, J)` as well as `I*J` to the quotient of `I` by `J`.

#### Tests on Ideals

The usual tests `f in J`, `is_subset(I, J)`, and `I == J` are available.

#### Saturation

If ``Rloc`` is the localization of a multivariate polynomial ring ``R`` at a multplicative subset ``U`` of ``R``,
then the ideal theory of ``Rloc`` is a simplified version of the ideal theory of ``R`` (see, for instance, [Eis95](@cite)).
In particular, each ideal ``I`` of ``Rloc`` is the extension $J\cdot Rloc$ of an ideal $J$ of $R$. The ideal

$$\{f\in R \mid uf\in J \text{ for some } u\in U\}$$

is independent of the choice of $J$ and is the largest ideal of ``R`` which extends to ``I``. It is, thus,
the contraction of ``I`` to ``R``,  that is, the preimage of ``I``  under the localization map.
We call this ideal the *saturation of ``I`` over ``R``*.  In OSCAR, it is obtained as follows:

```@docs
saturated_ideal(I::MPolyLocalizedIdeal)
```

## A Framework for Localizing Rings

For the convenience of the developer, we outline a general framework for creating concrete instances of localized rings in OSCAR,
addressing relevant abstract types as well as a standardized set of functions whose concrete behaviour must be implemented.

We roughly follow the outline of the previous subsection on localizing multivariate rings which provides illustrating examples.
With regard to notation, the name `Rloc` will refer to the localization of a commutative ring `R` with 1.

### Localized Rings

All multiplicatively closed subsets should belong to the `AbsMultSet` abstract type and all
localized rings should belong to the `AbsLocalizedRing` abstract type.

The basic functionality that has to be realized for any concrete instance of `AbsMultSet`
is the containment check for elements in multiplicatively closed subsets via the `in` function.

For each concrete instance of `AbsLocalizedRing`, the `Localization` constructor as well as the
functions `base_ring` and `inverted_set` need to be implemented. Moreover, as for any other type
of rings in OSCAR, methods for the standardized set of functions of OSCAR's
general [Ring Interface](@ref) must be supplied.

### Elements of Localized Rings

All elements of localized rings should belong to the `AbsLocalizedRingElem` abstract type.

Coercing (pairs of) elements of `R` into fractions in `Rloc` must be possible as indicated below:

```
   (Rloc::AbsLocalizedRing)(f::RingElem)
   (Rloc::AbsLocalizedRing)(f::RingElem, g::RingElem; check::Bool=true)
```
   
The first constructor maps the element `f` of `R` to the fraction `f//1` in `Rloc`.
The second constructor takes a pair `f, g` of elements of `R` to the fraction `f//g`
in `Rloc`. The default `check = true` stands for testing whether `g` is an admissible
denominator. As this test is often expensive, it may be convenient
to set `check = false`.

For any concrete instance of type `AbsLocalizedRingElem`, methods for the functions 
`parent`, `numerator`, and `denominator` must be provided. Moreover,
if a cancellation function for the type of fractions under consideration is
not yet available, such a function should be implemented and named
`reduce_fraction`.

### Homomorphisms From Localized Rings

The abstract type for homomorphisms from localized rings is `MPolyLocalizedRingHom`.
For each concrete instance, the functions `domain` and `codomain` as well as `restricted_map`
must be realized. Here, the latter function is meant to return the composition with the
localization map.

### Ideals in Localized Rings

All ideals in localized rings belong to the abstract type `AbsLocalizedIdeal`.
For a concrete instance, the constructors to be implemented are:

```
   ideal(W::AbsLocalizedRing, f::AbsLocalizedRingElem)
   ideal(W::AbsLocalizedRing, v::Vector{LocalizedRingElemType}) where {LocalizedRingElemType<:AbsLocalizedRingElem}
```

The usual getter functions  `base_ring`, `gens`, `ngens`, and `gen`   must be realized.

Moreover, a method for ideal membership via the `in` function is required.

