```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Groups of automorphisms

```@docs
automorphism_group(G::GAPGroup)
```

The following functions are available for automorphisms, some of them similar to the corresponding functions for homomorphisms of groups.
```@docs
is_invariant(f::GAPGroupElem{AutomorphismGroup{T}}, H::GAPGroup) where T <: GAPGroup
restrict_automorphism(f::GAPGroupElem{AutomorphismGroup{T}}, H::GAPGroup, A=automorphism_group(H)) where T <: GAPGroup
induced_automorphism(f::GAPGroupHomomorphism, mH::GAPGroupHomomorphism)
hom(x::GAPGroupElem{AutomorphismGroup{T}}) where T <: GAPGroup
```

## [Inner automorphisms](@id inner_automorphisms)

OSCAR provides the following functions to handle inner automorphisms of a group.
```@docs
inner_automorphism(g::GAPGroupElem)
is_inner_automorphism(f::GAPGroupHomomorphism)
inner_automorphism_group(A::AutomorphismGroup{T}) where T <: GAPGroup
```
