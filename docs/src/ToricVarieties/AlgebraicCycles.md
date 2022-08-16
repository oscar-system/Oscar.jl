```@meta
CurrentModule = Oscar
```

```@contents
Pages = ["AlgebraicCycles.md"]
```

# The Chow ring

Algebraic cycles are formal linear sum of irreducible
subvarieies over the integers. Perse, algebraic cycles
do not admit a well-defined intersection product.

To see this, think of intersecting a non-trivial algebraic
cycle `C` with itself. Of course, in set theory we can
intersect `C` with itself and the result is again `C`.
However, for a well-defined intersection theory, we would
ask that the self-intersection of `C` is an algebraic
cycle of strictly smaller dimension.

In theory, this is resolved by saying that the
self-intersection of `C` is given by intersecting `C` with
a distinct algebraic cycle `D` which is obtained by moving
`C` a little bit. The general phrase for this is to "move
`C` in general position".

This leads to a famous notion of equivalence among algebraic
cycles, the so-called rational equivalence. The set of
equivalence classes of algebraic cycles together with the
intersection product then furnishes the Chow ring of
the variety in question.

For complete and simplicial toric varieties, many things are
known about the Chow ring and algebraic cycles (cf. section 12.5
in [CLS11](@cite):
* By therorem 12.5.3 of [CLS11](@cite), there is an isomorphism
among the Chow ring, the cohomology ring and a certain quotient
of the non-graded Cox ring. (The Cox ring is perse graded by the
class group, however the ideal of linear relations that one has
to divide by in order to obtain the Chow ring is not necessarily
homogeneous with respect to this grading.)
* By lemma 12.5.1 of [CLS11](@cite), generators of the equivalence
classes of algebraic cycles are one-to-one to the cones in the fan
of the toric variety.


## Constructors

### General constructors

```@docs
RationalEquivalenceClass(v::AbstractNormalToricVariety, coefficients::Vector{T}) where {T <: IntegerUnion}
```

### Special constructors

```@docs
RationalEquivalenceClass(d::ToricDivisor)
RationalEquivalenceClass(c::ToricDivisorClass)
RationalEquivalenceClass(l::ToricLineBundle)
RationalEquivalenceClass(cc::CohomologyClass)
RationalEquivalenceClass(cc::ClosedSubvarietyOfToricVariety)
```

### Addition, subtraction and scalar multiplication

Algebraic cycles can be added and subtracted via the usual `+` and `-`
operators. Moreover, multiplication by scalars from the left is supported
for scalars which are integers or of type `fmpz`.

Note that one can easily define the Chow ring also a formal linear sums of
irreducible subvarieties with coefficients being rational numbers. We
support this more general ring and therefore also allow for left
multiplication with scalars of type `fmpq`.

### Intersection product

The intersection product of algebraic cycles is implemented via `*`.
This makes sense, since algebraic cycles on toric varieties are
elements of the Chow ring, which in turn is (a certain) quotient of
the Cox ring. Hence, internally, an algebraic cycle can be thought
of as a polynomial in this ring and the intersection product
corresponds to the product of two (equivalence classes of) polynomials.

An algebraic cycle can be intersected `n`- with itself via `^n`,
where `n` can be an integer of of type `fmpz`.

A closed subvarieties defines in a natural way a rational equivalence
class (cf. section on special constructors above). This allows to
compute intersection products among closed subvarieties and rational
equivalence classes in the Chow ring.


## Attributes

### Defining attributes

```@docs
toric_variety(ac::RationalEquivalenceClass)
polynomial(ac::RationalEquivalenceClass)
polynomial(ring::MPolyQuo, ac::RationalEquivalenceClass)
```

### Representant

In order to see a geometric interpretation of rational equivalence
classes of algebraic cycles most efficiently, it is best to replace
self-intersections by transverse complete intersections. Indeed,
within the regime of simplicial, complete toric varieties this is
always possible. However, this involves a choice. Consequently,
the following methods will pick a special choice and return
values for that particular choice of representant of the rational
equivalence class in question.

```@docs
representant(ac::RationalEquivalenceClass)
```
It can be rather convenient, to investigate such a representant in
order to understand the geometric meaning of a rational equivalence
class. For this purpose, we support the following methods.
```@docs
coefficients(ac::RationalEquivalenceClass)
components(ac::RationalEquivalenceClass)
```


### Other attributes

```@docs
cohomology_class(ac::RationalEquivalenceClass)
```


## Properties

One can check if a rational equivalence class of algebraic cycles
is trivial via `is_trivial`. Equality can be tested with `==`.


## Special attributes of toric varieties

```@docs
chow_ring(v::AbstractNormalToricVariety)
gens_of_rational_equivalence_classes(v::AbstractNormalToricVariety)
map_gens_of_chow_ring_to_cox_ring(v::AbstractNormalToricVariety)
```
