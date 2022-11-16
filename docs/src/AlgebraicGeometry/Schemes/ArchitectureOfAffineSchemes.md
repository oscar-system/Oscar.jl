```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["ArchitectureOfAffineSchemes.md"]
```


# Architecture of affine schemes

## Requirements on rings and ideals

Any type of ring ``R`` to be used within the schemes framework
must come with its own ideal type `IdealType<:Ideal` for which
we require the following interface to be implemented:
```
    # constructor of ideals in R
    ideal(R::RingType, g::Vector{<:RingElem})::Ideal

    # constructor for quotient rings
    quo(R::RingType, I::IdealType)::Tuple{<:Ring, <:Map}

    # ideal membership test
    in(f::RingElem, I::IdealType)::Bool

    # a (fixed) set of generators
    gens(I::IdealType)::Vector{<:RingElem}

    # writing an element as linear combination of the generators
    coordinates(f::RingElem, I::IdealType)::Vector{<:RingElem}
```
The latter function must return a vector ``v = (a_1,\dots, a_r)``
of elements in ``R`` such that ``f = a_1 \cdot g_1 + \dots + a_r \cdot g_r``
where ``g_1,\dots,g_r`` is the set of `gens(I)`. When ``f`` does
not belong to ``I``, it must error. Note that the ring returned by
`quo` must again be admissible for the `AbsSpec` interface.

With a view towards the use of the `ambient_coordinate_ring(X)` for computations,
it is customary to also implement
```
    saturated_ideal(I::IdealType)::MPolyIdeal
```
returning an ideal ``J`` in the `ambient_coordinate_ring(X)` with the property
that ``a \in I`` for some element ``a \in R`` if and only if
`lifted_numerator(a)` is in ``J``.


## Interplay between ambient coordinate ring and coordinate ring

Let ``X`` be an affine variety.
In practice, all computations in the coordinate ring `R = OO(X)` will be deferred to
computations in `P = ambient_coordinate_ring(X)` in one way or the other;
this is another reason to include the ambient affine space in our abstract
interface for affine schemes. In order to make the `ambient_coordinate_ring(X)`
accessible for this purpose, we need the following methods to be implemented
for elements ``a\in R`` of type `RingElemType`:
```
   lifted_numerator(a::RingElemType)
   lifted_denominator(a::RingElemType)
```
These must return representatives of the numerator and the denominator
of ``a``. Note that the denominator is equal to `one(P)` in case
``R \cong P`` or ``R \cong P/I``.

Recall that the coordinates ``x_i`` of ``X`` are induced by the coordinates of
the ambient affine space.
Moreover, we will assume that for homomorphisms from ``R``
there is a method
```
    hom(R::RingType, S::Ring, a::Vector{<:RingElem})
```
where `RingType` is the type of ``R`` and `a` the images
of the coordinates ``x_i`` in ``S``. This will be important
when we come to morphisms of affine schemes below.

Algebraically speaking, embedding the affine scheme ``X = Spec(R)`` over ``𝕜``
into an affine space with coordinate ring ``P = 𝕜[x₁,…,xₙ]`` corresponds to
a morphism of rings ``P → R`` with the property that for every other (commutative)
ring ``S`` and any homomorphism ``φ : R → S`` there is a morphism
``ψ : P → S`` factoring through ``φ`` and such that ``φ``
is uniquely determined by ``ψ``.
Note that the morphism ``P → R`` is induced by natural coercions.


## Existing types of affine schemes and how to derive new types

The abstract type for affine schemes is
```@docs
    AbsSpec{BaseRingType, RingType<:Ring}
```
For any concrete instance of this type, we require the following
functions to be implemented:
- `base_ring(X::AbsSpec)`,
- `OO(X::AbsSpec)`.

A concrete instance of this type is
```@docs
    Spec{BaseRingType, RingType}
```
It provides an implementation of affine schemes for rings ``R`` of type
`MPolyRing`, `MPolyQuo`, `MPolyLocalizedRing`, and `MPolyQuoLocalizedRing`
defined over the integers or algebraic field extensions of ``\mathbb Q``.
This minimal implementation can be used internally, when deriving new
concrete types `MySpec<:AbsSpec` such as, for instance,
group schemes, toric schemes, schemes of a particular dimension
like curves and surfaces, etc. To this end, one has to store
an instance `Y` of `Spec` in `MySpec` and implement the methods
```
    underlying_scheme(X::MySpec)::Spec # return Y as above
```
Then all methods implemented for `Spec` are automatically
forwarded to any instance of `MySpec`.

**Note:** The above method necessarily returns an instance of `Spec`!
Of course, it can be overwritten for any higher type `MySpec<:AbsSpec` as needed.


## Existing types of affine scheme morphisms and how to derive new types

Any abstract morphism of affine schemes is of the following type:
```@docs
    AbsSpecMor{DomainType<:AbsSpec,
               CodomainType<:AbsSpec,
               PullbackType<:Hecke.Map,
               MorphismType,
               BaseMorType
               }
```
Any such morphism has the attributes `domain`, `codomain` and `pullback`.
A concrete and minimalistic implementation exist for the type `SpecMor`:
```@docs
    SpecMor{DomainType<:AbsSpec,
            CodomainType<:AbsSpec,
            PullbackType<:Hecke.Map
           }
```
This basic functionality consists of
- `compose(f::AbsSpecMor, g::AbsSpecMor)`,
- `identity_map(X::AbsSpec)`,
- `restrict(f::AbsSpecMor, X::AbsSpec, Y::AbsSpec; check::Bool=true)`,
- `==(f::AbsSpecMor, g::AbsSpecMor)`,
- `preimage(f::AbsSpecMor, Z::AbsSpec)`.
In particular, for every concrete instance of a type `MySpec<:AbsSpec` that
implements `underlying_scheme`, this basic functionality of `SpecMor`
should run naturally.

We may derive higher types of morphisms of affine schemes `MySpecMor<:AbsSpecMor`
by storing an instance `g` of `SpecMor` inside an instance `f` of
`MySpecMor` and implementing
```
    underlying_morphism(f::MySpecMor)::SpecMor # return g
```
For example, this allows us to define closed embeddings.
