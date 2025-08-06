```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Morphisms of projective schemes

Let ``Q = B[y_0, \dots, y_n]/J`` and ``P = A[x_0,\dots,x_m]/I`` be 
graded affine algebras over `base_ring`s `A` and `B`, respectively. 
A morphism ``\varphi : \mathrm{Proj}(Q) \to \mathrm{Proj}(P)`` is modeled 
via a morphism of graded algebras ``\varphi^* : P \to Q``. 
In the case of `A != B`, this involves a non-trivial morphism 
of rings ``A \to B``.

## Abstract types and basic interface 
At the moment we have no abstract type for such morphisms and no interface spelled 
out. 

## Types 
```@docs
ProjectiveSchemeMor
```

## Constructors
```@docs
morphism(P::AbsProjectiveScheme, Q::AbsProjectiveScheme, f::Map; check::Bool=true )
morphism(P::AbsProjectiveScheme, Q::AbsProjectiveScheme, f::Map, h::SchemeMor; check::Bool=true )
morphism(X::AbsProjectiveScheme, Y::AbsProjectiveScheme, a::Vector{<:RingElem})
```
## Attributes
As every instance of `Map`, a morphism of projective schemes can be asked for its (co-)domain:
```julia
domain(phi::ProjectiveSchemeMor) 
codomain(phi::ProjectiveSchemeMor)
```
Moreover, we provide getters for the associated morphisms of rings:
```@docs
pullback(phi::ProjectiveSchemeMor)
base_ring_morphism(phi::ProjectiveSchemeMor) 
base_map(phi::ProjectiveSchemeMor)
map_on_affine_cones(phi::ProjectiveSchemeMor)
```
## Methods
```@docs
covered_scheme_morphism(f::AbsProjectiveSchemeMorphism)
```

