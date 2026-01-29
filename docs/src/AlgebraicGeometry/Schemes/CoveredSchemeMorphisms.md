```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Morphisms of covered schemes

Suppose ``f : X \to Y`` is a morphism of `AbsCoveredScheme`s. Theoretically, and hence 
also technically, the required information behind ``f`` is a list of morphisms 
of affine schemes ``f_i : U_i \to V_{F(i)}`` for some pair of `Covering`s ``\left\{U_i\right\}_{i \in I}``
of ``X`` and ``\left\{V_j\right\}_{j \in J}`` of ``Y`` and a map of indices ``F : I \to J``
This information is held by a `CoveringMorphism`:
```@docs 
    CoveringMorphism
```
The basic functionality of `CoveringMorphism`s comprises `domain` and `codomain` which 
both return a `Covering`, together with 
```julia
getindex(f::CoveringMorphism, U::AbsAffineScheme)
```
which for ``U = U_i`` returns the `AbsAffineSchemeMor` ``f_i : U_i \to V_{F(i)}``.

Note that, in general, neither the `domain` nor the `codomain` of the `covering_morphism` of 
`f : X \to Y` need to coincide with the `default_covering` of ``X``, respectively ``Y``. 
In fact, one will usually need to restrict to a refinement of the `default_covering` of ``X``
in order to realize the covering morphism in the first place. 

## The interface for morphisms of covered schemes
Every `AbsCoveredSchemeMorphism` ``f : X \to Y`` is required to implement the following minimal 
interface.
```julia
domain(f::AbsCoveredSchemeMorphism)                 # returns X
codomain(f::AbsCoveredSchemeMorphism)               # returns Y
covering_morphism(f::AbsCoveredSchemeMorphism)      # returns the underlying covering morphism {f_i}
```
For the user's convenience, also the domain and codomain 
of the underlying `covering_morphism` are forwarded as `domain_covering` and 
`codomain_covering`, respectively, together with `getindex(phi::CoveringMorphism, U::AbsAffineScheme)` 
as `getindex(f::AbsCoveredSchemeMorphism, U::AbsAffineScheme)`.
    
The minimal concrete type of an `AbsCoveredSchemeMorphism` which 
implements this interface, is `CoveredSchemeMorphism`.

## Special types of morphisms of covered schemes
### Closed embeddings
```@docs 
CoveredClosedEmbedding
image_ideal(phi::CoveredClosedEmbedding)
```

### Composite morphisms
```@docs
CompositeCoveredSchemeMorphism
```

### Morphisms from rational functions

Suppose ``X`` and ``Y`` are two irreducible and reduced varieties. Then a morphism 
``f : X \to Y`` might be given by means of the following data. Let ``V \subset Y`` be 
some dense affine patch with `coordinates` ``x_1,\dots,x_n`` (i.e. the `gens` of `OO(V)`). 

These ``x_i`` extend to rational functions ``v_i`` on ``Y`` and these pull back to rational functions 
``f^* v_i = u_i`` on ``X``. On every affine patch ``U \subset X`` there now exists 
some maximal Zariski-open subset ``W \subset U`` (which need not be affine), such that 
all the ``f^* v_i`` extend to regular functions on ``W``. Hence, one can realize the 
morphisms of affine schemes ``f_j : W_j \to V`` for some open covering of ``W``. 

Similarly, for every other non-empty patch ``V_2`` of ``Y`` the pullback of `gens(OO(V_2))` can 
be computed from the ``f^* v_i`` and extended maximally to some ``W \subset U`` for every patch 
``U`` of ``X``. Altogether, this allows to compute a full `CoveringMorphism` and 
-- at least in theory -- an instance of `CoveredSchemeMorphism`. In practice, 
however, this computation is usually much too expensive to really be carried out, while 
the data necessary to compute various pullbacks and/or pushforwards of objects defined 
on ``X`` and ``Y`` can be extracted from the ``f^* v_i`` more directly. 

A lazy concrete data structure to house this kind of morphism is 
```@docs
MorphismFromRationalFunctions
```
Note that the key idea of this data type is to *not* use the `underlying_morphism` 
together with its `covering_morphism`, but to find cheaper ways to do computations! 
The computation of the `underlying_morphism` is triggered by any call to functions 
which have not been overwritten with a special method for `f::MorphismFromRationalFunctions`. 
However, this computation should be considered as way too expensive besides some small examples. 

For instance, if one wants to pull back a prime ideal sheaf ``\mathcal I`` on ``Y`` along 
some *isomorphism* ``f : X \to Y``, then one only needs to find one realization 
``f_j : U_j \to V_{F(j)}`` of ``f`` on affine patches ``U_j`` of ``X`` and ``V_{F(j)}`` 
of ``Y`` such that ``\mathcal I(V_{F(j)})\neq 0`` and ``f_j^* \mathcal I(V_{F(j)}) \neq 0``. 
Then ``f^* \mathcal I`` can be extended uniquely to all of ``X`` from ``U_j`` and 
there is no need to realize the full `covering_morphism` of ``f``.
In order to facilitate such computations as lazy as possible, there are various fine-grained 
entry points and caching mechanisms to realize ``f`` on open subsets:
```@docs
realize_on_patch(Phi::MorphismFromRationalFunctions, U::AbsAffineScheme)
realize_on_open_subset(Phi::MorphismFromRationalFunctions, U::AbsAffineScheme, V::AbsAffineScheme)
realization_preview(Phi::MorphismFromRationalFunctions, U::AbsAffineScheme, V::AbsAffineScheme)
random_realization(Phi::MorphismFromRationalFunctions, U::AbsAffineScheme, V::AbsAffineScheme)
cheap_realization(Phi::MorphismFromRationalFunctions, U::AbsAffineScheme, V::AbsAffineScheme)
realize_maximally_on_open_subset(Phi::MorphismFromRationalFunctions, U::AbsAffineScheme, V::AbsAffineScheme)
realize(Phi::MorphismFromRationalFunctions)
```

