```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

# General schemes and their interfaces

Arbitrary schemes over a base ring ``\mathbb k`` which are given by means 
of their affine patches and glueings are instances of the abstract type
```@docs
Scheme{BaseRingType<:Ring}
```
Morphisms of schemes shall be derived from the abstract type
```@docs
    SchemeMor{DomainType, CodomainType, MorphismType, BaseMorType}
```

## Abstract affine schemes

### The mathematical interface
Let ``\mathbb k`` be a commutative noetherian base ring 
(in practice: an algebraic extension of ``\mathbb Q`` or ``\mathbb F_p``). 
An affine scheme ``X = \mathrm{Spec}(R)`` over ``\mathbb k`` is an instance of 
the abstract type
```@docs 
    AbsSpec{BaseRingType, RingType<:Ring}
```
For any concrete instance of this type, we require the following interface 
to be implemented. 
```@docs
    base_ring(X::AbsSpec) 
    OO(X::AbsSpec)
```

### The algorithmic interface
For most affine schemes ``X = \mathrm{Spec}(R)`` 
over ``\mathbb k``, there is a 'governing' polynomial 
``\mathbb k``-algebra ``P = \mathbb{k}[x_1,\dots,x_n]`` 
in the following sense:
```@docs
    ambient_ring(X::AbsSpec)
```
For instance, this is the case whenever ``R`` is a quotient 
ring of ``P``, a localization of ``P``, or a localization 
of a quotient ring of ``P``; but also for 
power series rings in multiple variables. 

In practice, 
all computations in `OO(X)` will be deferred to computations in 
`ambient_ring(X)` in one way or the other; that is another 
reason to include this getter in our abstract interface for 
affine schemes. In order to make the `ambient_ring(X)` accessible 
for this purpose, we need the following methods to be implemented 
for elements ``a\in R`` of type `RingElemType`:
```
   lifted_numerator(a::RingElemType) 
   lifted_denominator(a::RingElemType)
```
These must return representatives of the numerator and the denominator 
of ``a``. Note that the denominator is equal to `one(P)` in case 
``R \cong P`` or ``R \cong P/I``. 

Given the above characterization of the `ambient_ring(X)`,
we shall refer to the ``x_i`` above as the *coordinates* of ``X``.
Moreover, we will assume that for homomorphisms from ``R`` 
there is a a method 
```
    hom(R::RingType, S::Ring, a::Vector{<:RingElem})
```
where `RingType` is the type of ``R`` and `a` the images 
of the coordinates ``x_i`` in ``S``. This will be important 
when we come to morphisms of affine schemes below.

### Required interface for the ideals to be used
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

With a view towards the use of the `ambient_ring(X)` for computations, 
it is customary to also implement 
```
    saturated_ideal(I::IdealType)::MPolyIdeal
```
returning an ideal ``J`` in the `ambient_ring(X)` with the property 
that ``a \in I`` for some element ``a \in R`` if and only if 
`lifted_numerator(a)` is in ``J``.

## A minimal implementation of the `AbsSpec` interface
An implementation that currently supports rings ``R`` of type 
`MPolyRing`, `MPolyQuo`, `MPolyLocalizedRing`, and `MPolyQuoLocalizedRing`
defined over the integers or algebraic field extensions of ``\mathbb Q``
is given by
```@docs
    Spec{BaseRingType, RingType}
```
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
In particular, the following additional package of 
commands will run out of the box:
```@docs
  subscheme(X::AbsSpec, f::RingElem)
  hypersurface_complement(X::AbsSpec, f::RingElem)
```
Schemes can be compared based on their `ambient_ring`, leading to the methods for 
```@docs 
    issubset(X::AbsSpec, Y::AbsSpec)
    is_open_embedding(X::AbsSpec, Y::AbsSpec)
    is_closed_embedding(X::AbsSpec, Y::AbsSpec)
```
Moreover, we have
```@docs 
    closure(X::AbsSpec, Y::AbsSpec) 
    product(X::AbsSpec, Y::AbsSpec)
```
**Note:** The above methods necessarily return instances of `Spec`!
Of course, they can again be overwritten for any higher 
type `MySpec<:AbsSpec` as needed. 


## Morphisms of affine schemes

### Abstract interface for morphisms of affine schemes
To this end we have the abstract type
```@docs
    AbsSpecMor{DomainType<:AbsSpec, 
               CodomainType<:AbsSpec, 
               PullbackType<:Hecke.Map,
               MorphismType, 
               BaseMorType
               }
```
The interface for this abstract type is 
```@docs
    domain(f::AbsSpecMor)
    codomain(f::AbsSpecMor)
    pullback(f::AbsSpecMor)
```

### A minimal implementation of the `AbsSpecMor` interface
Again, there is a minimal implementation of this interface for all instances 
of `AbsSpec` that meet the above requirements on the involved 
rings (especially the required methods for `hom`):
```@docs 
    SpecMor{DomainType<:AbsSpec, 
            CodomainType<:AbsSpec, 
            PullbackType<:Hecke.Map
           }
```
The default constructor is
```
    SpecMor(X::AbsSpec, Y::AbsSpec, a::Vector{<:RingElem}; check::Bool=true)
```
where one has to specify the domain and codomain of the map, and 
the images `a` of the coordinates (the generators of `ambient_ring(Y)`) 
under the pullback map ``\mathscr{O}(Y) \to \mathscr{O}(X)``. Expensive checks can be 
turned off by setting `check=false`.

The basic functionality for `SpecMor` comprises
```
   # the composition 
   compose(f::AbsSpecMor, g::AbsSpecMor) 

   # the identity map
   identity_map(X::AbsSpec)

   # restriction to subsets in the domain and codomain
   restrict(f::AbsSpecMor, X::AbsSpec, Y::AbsSpec; check::Bool=true)

   # equality test
   ==(f::AbsSpecMor, g::AbsSpecMor)

   # preimages of subschemes in the codomain
   preimage(f::AbsSpecMor, Z::AbsSpec)
```
In particular, for every concrete instance of a type `MySpec<:AbsSpec` that 
implements `underlying_scheme`, this basic functionality of `SpecMor` 
should run naturally.

### Some particular higher types of morphisms

We may derive higher types of morphisms of affine schemes `MySpecMor<:AbsSpecMor` 
by storing an instance `g` of `SpecMor` inside an instance `f` of 
`MySpecMor` and implementing 
```
    underlying_morphism(f::MySpecMor)::SpecMor # return g
```
For example, this allows us to define closed embeddings:
```@docs
    ClosedEmbedding{DomainType, CodomainType, PullbackType}
```
In addition to the standard getters and methods for instances 
of `SpecMor`, we also have 
```@docs
    image_ideal(f::ClosedEmbedding)
```







