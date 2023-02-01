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

We recall the definition of localization. All rings considered are commutative,  with multiplicative identity 1.
Let ``R`` be a ring, and let ``U \subset R`` be a *multiplicatively closed subset.* That is,
```math
1 \in U  \;\text{ and }\;  u, v \in U \;\Rightarrow \; u\cdot v \in U.
```
Consider the equivalence relation on ``R\times U`` defined by setting
```math
(r,u)\sim (r', u') \;\text{ iff }\; v(r u'-u r')=0 \;{\text{ for some }}\; v\in U.
```
Write ``\frac{r}{u}`` for the equivalence class of ``(r, u)`` and ``R[U^{-1}]`` for the set of all equivalence classes.
Mimicking the standard arithmetic for fractions, ``R[U^{-1}]`` can be made into a ring. This ring is called the
*localization of ``R`` at ``U``*. It comes equipped with  the natural ring homomorphism
```math
\iota : R\rightarrow R[U^{-1}],\; r \mapsto \frac{r}{1}.
```
Given an ``R``-module ``M``, the analogous construction yields an ``R[U^{-1}]``-module ``M[U^{-1}]`` which is
called the *localization  of ``M`` at ``U``*.  See, for instance, [Eis95](@cite).

OSCAR supports several types of multiplicatively closed subsets of multivariate polynomial rings.
The related functionality is not only used to construct localizations of
- multivariate polynomial rings (OSCAR type `MPolyRing`),
but also to construct localizations of
- quotients of multivariate polynomial rings (OSCAR type `MPolyQuo`).

For the latter, recall that if ``A`` is a quotient of a multivariate polynomial ring ``R``,
``p: R \rightarrow A`` is the projection map, and ``S`` is a multiplicatively closed subset of ``A``,
then the preimage $U = p^{-1}(S)$ is a multiplicatively closed subset of ``R`` such that ``S = p(U)``.
Creating the localization of the ring ``A`` at ``S `` in OSCAR amounts to localize the ``R``-module
``A`` at ``U``, and equip the result with the appropriate ring structure. Hence, the corresponding
constructor takes as input ``A`` and ``U`` rather than ``A`` and ``S``.

!!! note
    Most functions described in this section rely on the computation of standard bases. Recall that OSCAR supports standard
    bases for multivariate polynomial rings over fields (exact fields supported by OSCAR) and for multivariate
	polynomial rings over the integers.	

## Types

The OSCAR types discussed in this section are all parametrized. To simplify the presentation,
details on the parameters are omitted.

All types for multiplicatively closed subsets of rings belong to the abstract type `AbsMultSet`.
For multiplicatively closed subsets of multivariate polynomial rings, there are the concrete descendants
`MPolyComplementOfKPointIdeal`, `MPolyComplementOfPrimeIdeal`, and `MPolyPowersOfElement`.

The general abstract type for localizations of rings is `AbsLocalizedRing`. For localizations of multivariate
polynomial rings, there is the concrete subtype `MPolyLocalizedRing`. For localizations of quotients of
multivariate polynomial rings, there is the concrete subtype `MPolyQuoLocalizedRing`. 


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

Containment in multiplicatively closed subsets can be checked via the `in` function:

```@docs
in(f::MPolyElem, U::AbsMPolyMultSet)
```

### Localized Rings

```@docs
localization(R::MPolyRing, U::AbsMPolyMultSet)
```

```@docs
localization(RQ::MPolyQuo, U::AbsMPolyMultSet)
```
	
## Data associated to Localized Rings

If `Rloc` is the localization of a multivariate polynomial ring `R`  at a multiplicatively closed subset
`U` of `R`, then

- `base_ring(Rloc)` refers to `R`, and
- `inverted_set(Rloc)` to `U`.

If `RQ` is a quotient of a multivariate polynomial ring `R`, `p : R -> RQ` is the projection map, `U`  is a
multiplicatively closed subset of `R`, and `RQL` is the localization of `RQ` at `p(U)`, then

- `base_ring(RQL)` refers to `R`, and
- `inverted_set(RQL)` to `U`.

This reflects the internal way of creating localizations of quotients of multivariate polynomial rings.

##### Examples

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> P = ideal(R, [x])
ideal(x)

julia> U = complement_of_ideal(P)
complement of ideal(x)

julia> Rloc, _ = localization(U);

julia> R === base_ring(Rloc)
true

julia> U === inverted_set(Rloc)
true
```

```jldoctest
julia> T, t = PolynomialRing(QQ, "t");

julia> K, a =  NumberField(2*t^2-1, "a");

julia> R, (x, y) = PolynomialRing(K, ["x", "y"]);

julia> I = ideal(R, [2*x^2-y^3, 2*x^2-y^5])

ideal(2*x^2 - y^3, 2*x^2 - y^5)

julia> P = ideal(R, [y-1, x-a])
ideal(y - 1, x - a)

julia> U = complement_of_ideal(P)
complement of ideal(y - 1, x - a)

julia> RQ, _ = quo(R, I);

julia> RQL, _ = localization(RQ, U);

julia> R == base_ring(RQL)
true

julia> U == inverted_set(RQL)
true
```

## Elements of Localized Rings

### Types

The general abstract type for elements of localizations of rings is `AbsLocalizedRingElem`.
For elements of localizations of multivariate polynomial rings, there is the concrete subtype `MPolyLocalizedRingElem`.
For elements of localizations of quotients of multivariate polynomial rings, there is the concrete subtype `MPolyQuoLocalizedRingElem`.

### Creating Elements of Localized Rings

If `Rloc` is the localization of a multivariate polynomial ring `R`  at a multiplicatively closed subset
`U` of `R`, then elements of `Rloc` are created as (fractions of) images of elements of `R` under
the localization map or by coercing (pairs of) elements of `R` into fractions. 

If `RQ` is a quotient of a multivariate polynomial ring `R`, `p : R -> RQ` is the projection map, `U`  is a
multiplicatively closed subset of `R`, and `RQL` is the localization of `RQ` at `p(U)`, then elements of
`RQL` are created similarly, starting from elements of `R`.

##### Examples

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> P = ideal(R, [x])
ideal(x)

julia> U = complement_of_ideal(P)
complement of ideal(x)

julia> Rloc, iota = localization(U);

julia> iota(x)
x//1

julia> Rloc(x)
x//1

julia> f = iota(y)//iota(z)
y//z

julia> g = Rloc(y, z)
y//z

julia> X, Y, Z = Rloc.(gens(R));

julia> h = Y//Z
y//z

julia> f == g == h
true

julia> f+g
(2*y)//z

julia> f*g
y^2//z^2
```

```jldoctest
julia> T, t = PolynomialRing(QQ, "t");

julia> K, a =  NumberField(2*t^2-1, "a");

julia> R, (x, y) = PolynomialRing(K, ["x", "y"]);

julia> I = ideal(R, [2*x^2-y^3, 2*x^2-y^5])
ideal(2*x^2 - y^3, 2*x^2 - y^5)

julia> P = ideal(R, [y-1, x-a])
ideal(y - 1, x - a)

julia> U = complement_of_ideal(P)

complement of ideal(y - 1, x - a)

julia> RQ, p = quo(R, I);

julia> RQL, iota = Localization(RQ, U);

julia> phi = compose(p, iota);

julia> phi(x)
x//1

julia> RQL(x)
x//1

julia> f = phi(x)*inv(phi(y))
x//y

julia> g = RQL(x, y)
x//y

julia> X, Y = gens(RQL);

julia> h = X*inv(Y)
x//y

julia> f == g == h
true

julia> f+g
(2*x)//y

julia> f*g
x^2//y^2
```

### Data Associated to Elements of Localized Rings

If `Rloc` is a localization of a multivariate polynomial ring `R`, and `f` is an element of `Rloc`, internally
represented by a pair `(r, u)` of elements of `R`, then 
- `parent(f)` refers to `Rloc`,
- `numerator(f)` to `r`, and
- `denominator(f)` to `u`.
If `RQL` is a localization of a quotient `RQ` of a multivariate polynomial ring `R`, and `f` is an element of `RQL`,
internally represented by a pair `(r, u)` of elements of `R`, then
- `parent(f)` refers to `Rloc`,
- `numerator(f)` to the image of `r` in `RQ`, and
- `denominator(f)` to the image of `u` in `RQ`.
That is, the behaviour of the functions `numerator` and `denominator` reflects our mathematical viewpoint
of representing `f` by pairs of elements of `RQ`, rather than reflecting the internal representation of `f`.

##### Examples

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> P = ideal(R, [x])
ideal(x)

julia> U = complement_of_ideal(P)
complement of ideal(x)

julia> Rloc, iota = localization(U);

julia> f = iota(x)//iota(y)
x//y

julia> parent(f)
localization of Multivariate Polynomial Ring in x, y, z over Rational Field at the complement of ideal(x)

julia> g = iota(y)//iota(z)
y//z

julia> r = numerator(f*g)
x

julia> u = denominator(f*g)
z

julia> typeof(r) == typeof(u) <: MPolyElem
true
```

```jldoctest
julia> T, t = PolynomialRing(QQ, "t");

julia> K, a =  NumberField(2*t^2-1, "a");

julia> R, (x, y) = PolynomialRing(K, ["x", "y"]);

julia> I = ideal(R, [2*x^2-y^3, 2*x^2-y^5])
ideal(2*x^2 - y^3, 2*x^2 - y^5)

julia> P = ideal(R, [y-1, x-a])
ideal(y - 1, x - a)

julia> U = complement_of_ideal(P)
complement of ideal(y - 1, x - a)

julia> RQ, p = quo(R, I);

julia> RQL, iota = Localization(RQ, U);

julia> phi = compose(p, iota);

julia> f = phi(x)
x//1

julia> parent(f)
Localization of Quotient of Multivariate Polynomial Ring in x, y over K by ideal(2*x^2 - y^3, 2*x^2 - y^5) at the multiplicative set complement of ideal(y - 1, x - a)

julia> g = f*inv(phi(y))
x//y

julia> r = numerator(f*g)
x^2

julia> u = denominator(f*g)
y

julia> typeof(r) == typeof(u) <: MPolyQuoElem
true
```

### Tests on Elements of Localized Rings

```@docs
is_unit(f::MPolyLocalizedRingElem)
```

```@docs
is_unit(f::MPolyQuoLocalizedRingElem)
```

## Homomorphisms from Localized Rings

The general abstract type for ring homomorphisms starting from localized rings is `AbsLocalizedRingHom`.
For ring homomorphisms starting from localizations of multivariate polynomial rings, there is the concrete
subtype `MPolyLocalizedRingHom`. For ring homomorphisms starting from quotients of multivariate polynomial
rings, there is the concrete subtype `MPolyQuoLocalizedRingHom`. We describe the construction of such
homomorphisms. Let
- `R` be a multivariate polynomial ring
- `U` be a multiplicatively closed subset  of `R`,
- `RQ` be a quotient of `R` with projection map `p: R -> RQ`,
- `Rloc` (`RQL`) be the localization of `R` at `U` (of `RQ` at `p(U)`), and
- `S` be another ring.
Then, to give a ring homomorphism `PHI`  from `Rloc` to `S` (from`RQL` to `S`) is the same
as to give a ring homomorphism `phi` from `R` to `S` which sends elements of `U` to units
in `S`. That is, `PHI` is determined by the composition of `PHI` with the localization map
`R -> Rloc` (with the composition `R -> RQ -> RQL` of the localization map with the
projection map). The constructors below take this into account.

```@docs
hom(Rloc::MPolyLocalizedRing, S::Ring, F::Map)
```
	
Given a ring homomorphism `PHI` from `Rloc` to `S` (from `RQL` to `S`), `domain(PHI)` and `codomain(PHI)`
refer to `Rloc` and `S` (`RQL`  and `S`), respectively. The corresponding homomorphism `phi` from `R`
to `S` is recovered as follows:

```@docs
restricted_map(PHI::MPolyLocalizedRingHom)
```

## Ideals in Localized Rings
	  
### Types

The general abstract type for ideals in localized rings is `AbsLocalizedIdeal`. For ideals in  localizations of multivariate polynomial rings,
there is the concrete subtype `MPolyLocalizedIdeal`. For ideals in localizations of quotients of multivariate polynomial rings, there is
the concrete subtype `MPolyQuoLocalizedIdeal`.

### Constructors

Given a localization `Rloc` of a multivariate polynomial ring `R`, and given a vector `V` of elements of
`Rloc` (of `R`),  the ideal of `Rloc` which is generated by (the images) of the entries of `V`
is created by entering `ideal(Rloc, V)`. Similarly for ideals in localizations of quotients of
multivariate polynomial rings.

##### Examples

```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"]);

julia> f = x^3+y^4
x^3 + y^4

julia> V = [derivative(f, i) for i=1:2]
2-element Vector{fmpq_mpoly}:
 3*x^2
 4*y^3

julia> U = complement_of_ideal(R, [0, 0]);

julia> Rloc, _ = localization(R, U);

julia> MI = ideal(Rloc, V)
ideal in localization of Multivariate Polynomial Ring in x, y over Rational Field at the complement of maximal ideal corresponding to point with coordinates fmpq[0, 0] generated by the elements (3*x^2)//1, (4*y^3)//1
```

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

If ``Rloc`` is the localization of a multivariate polynomial ring ``R`` at a multiplicative subset ``U`` of ``R``,
then the ideal theory of ``Rloc`` is a simplified version of the ideal theory of ``R`` (see, for instance, [Eis95](@cite)).
In particular, each ideal ``I`` of ``Rloc`` is the extension $J\cdot Rloc$ of an ideal $J$ of $R$. The ideal

$$\{f\in R \mid uf\in J \text{ for some } u\in U\}$$

is independent of the choice of $J$ and is the largest ideal of ``R`` which extends to ``I``. It is, thus,
the contraction of ``I`` to ``R``,  that is, the preimage of ``I``  under the localization map.
We call this ideal the *saturation of ``I`` over ``R``*.  In OSCAR, it is obtained by entering
`saturated_ideal(I)`.

If ``RQL`` is the localization of a quotient ``RQ`` of a multivariate polynomial ring ``R``, and
``I`` is an ideal of ``RQL``, then the return value of `saturated_ideal(I)` is the preimage of
the saturation of ``I`` over ``RQ`` under the projection map ``R \rightarrow RQ`` (and not
the saturation of ``I`` over ``RQ`` itself).

```@docs
saturated_ideal(I::MPolyLocalizedIdeal)
```

