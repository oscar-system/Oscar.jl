```@meta
CurrentModule = Oscar 
DocTestSetup = quote
  using Oscar
end
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["localizations.md"]
```

# Localized Rings and Their Ideals

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

Most functions for handling localizations of multivariate polynomial rings rely on Gröbner bases. More precisely,
if ``R`` is a multivariate polynomial ring, Gröbner bases allow one to solve  ideal- and module-theoretic
questions concerning ``R[U^{-1}]`` by computations over ``R``.

!!! note
    Recall that in OSCAR, Gröbner bases are implemented for multivariate polynomial rings over fields (exact fields supported by OSCAR) and for multivariate polynomial rings over the integers.

!!! note
    In the local context, following Hironaka and Grauert, it is often common to use the name *standard basis* instead of Gröbner basis. See, for example, [GP08](@cite).
	

## Types

All OSCAR types for multiplicatively closed subsets of commutative rings with unit
belong to the abstract type `AbsMultSet`.  For multiplicatively
closed subsets of multivariate polynomial rings, there are the concrete descendants
`MPolyComplementOfKPointIdeal`, `MPolyComplementOfPrimeIdeal`, and `MPolyPowersOfElement`.

The general abstract type for localizations of commutative rings with unit is `AbsLocalizedRing`.
For localizations of multivariate polynomial rings, there is the concrete subtype `MPolyLocalizedRing`.

## Constructors

### Multiplicatively Closed Subsets

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

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> S = complement_of_ideal(R, [0, 0 ,0])
complement of maximal ideal corresponding to point with coordinates fmpq[0, 0, 0]

julia> y in S
false

julia> P = ideal(R, [x])
ideal(x)

julia> T = complement_of_ideal(P)
complement of ideal(x)

julia> y in T
true

julia> f = x
x

julia> U = powers_of_element(f)
powers of fmpq_mpoly[x]

julia> x^3 in U
true

julia> (1+y)*x^2 in product(S, U)
true

```


### Localized Rings

```@docs
Localization(U::AbsMPolyMultSet)
```
	
## Data associated to Localized Rings

If `Rloc` is the localization of a multivariate polynomial ring `R`  at a multiplicatively closed subset
`U` of `R`, then

- `base_ring(Rloc)` refers to `R`, and
- `inverted_set(Rloc)` to `U`.

###### Examples

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> P = ideal(R, [x])
ideal(x)

julia> U = complement_of_ideal(P)
complement of ideal(x)

julia> Rloc, _ = Localization(U);

julia> R === base_ring(Rloc)
true

julia> U === inverted_set(Rloc)
true

```

## Elements of Localized Rings

### Types

The general abstract type for elements of localizations of commutative rings with unit is `AbsLocalizedRingElem`.
For elements of localizations of multivariate polynomial rings, there is the concrete subtype `MPolyLocalizedRingElem`.


### Creating Elements of Localized Rings

Elements of a localized multivariate polynomial ring ``R[U^{-1}]`` are created as (fractions of) images of elements of $R$
under the localization map or by directly coercing (pairs of) elements of $R$ into fractions.

##### Examples

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> P = ideal(R, [x])
ideal(x)

julia> U = complement_of_ideal(P)
complement of ideal(x)

julia> Rloc, iota = Localization(U);

julia> f = iota(x)
x//1

julia> f == Rloc(x)
true

julia> g = Rloc(y, z)
y//z

julia> f+g
(x*z + y)//z

julia> f*g
(x*y)//z

```

### Data Associated to Elements of Localized Rings

Given an element `f` of a localized multivariate ring polynomial `Rloc`, 
- `parent(f)` refers to `Rloc`,
- `numerator(f)` to the numerator  of the internal representative of `f`, and
- `denominator(f)` to the denominator of the internal representative of `f`.

##### Examples

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> P = ideal(R, [x])
ideal(x)

julia> U = complement_of_ideal(P)
complement of ideal(x)

julia> Rloc, iota = Localization(U);

julia> f = iota(x)//iota(y)
x//y

julia> parent(f)
localization of Multivariate Polynomial Ring in x, y, z over Rational Field at the complement of ideal(x)

julia> g = iota(y)//iota(z)
y//z

julia> numerator(f*g)
x

julia> denominator(f*g)
z

```

## Homomorphisms from Localized Rings

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
	  
## Ideals in Localized Rings
	  
### Types

Ideals in localized rings are of concrete type `MPolyLocalizedIdeal`.

### Constructors

Given a localization `Rloc` of a multivariate polynomial ring `R`, and given a vector `V` of elements of
`Rloc` (of `R`),  the ideal of `Rloc` which is generated by (the images) of the entries of `V`
is created by entering `ideal(Rloc, V)`. 

### Data Associated to Ideals

If `I` is an ideal of a localized multivariate polynomial ring  `Rloc`, then

- `base_ring(I)` refers to `Rloc`,
- `gens(I)` to the generators of `I`,
- `ngens(I)` to the number of these generators, and
- `gen(I, k)` as well as `I[k]` to the `k`-th such generator.

### Operations on Ideals

If `I`,  `J` are ideals of a localized multivariate polynomial ring  `Rloc`, then

- `I^k` refers to the `k`-th  power of `I`, 
- `I+J`, `I*J`,  and `intersect(I, J)` to the sum, product, and intersection of `I` and  `J`, and
- `quotient(I, J)` as well as `I:J` to the ideal quotient of `I` by `J`.

### Tests on Ideals

The usual tests `f in J`, `issubset(I, J)`, and `I == J` are available.

### Saturation

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

