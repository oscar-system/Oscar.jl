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
Pages = ["subgroups.md"]
```

# [Subgroups](@id subgroups)

The following functions are available in OSCAR for subgroup properties:

```@docs
sub(G::GAPGroup, gens::AbstractVector{<:GAPGroupElem}; check::Bool = true)
is_subgroup
embedding(G::T, H::T) where T <: GAPGroup
index(G::T, H::T) where T <: GAPGroup
is_normal(G::T, H::T) where T <: GAPGroup
is_characteristic(G::T, H::T) where T <: GAPGroup
```

## Standard subgroups

The following functions are available in OSCAR to obtain standard subgroups of
a group `G`. Every such function returns a tuple `(H,f)`, where `H` is a group
of the same type of `G` and `f` is the embedding homomorphism of `H` into `G`.

```@docs
trivial_subgroup
center(G::GAPGroup)
sylow_subgroup(G::GAPGroup, p::IntegerUnion)
derived_subgroup
fitting_subgroup
frattini_subgroup
radical_subgroup
socle
pcore(G::GAPGroup, p::IntegerUnion)
intersect(V::T...) where T<:GAPGroup
```

The following functions return a vector of subgroups.

```@docs
subgroups(G::GAPGroup)
normal_subgroups
maximal_subgroups
maximal_normal_subgroups
minimal_normal_subgroups
characteristic_subgroups
derived_series
sylow_system
hall_subgroup_reps
hall_system
complement_system
```

!!! note
    When a function returns a vector of subgroups,
    the output consists in the subgroups only;
    the embeddings are not returned as well.
    To get the embedding homomorphism of the subgroup `H` in `G`,
    one can type `embedding(G,H)`.


## Conjugation action of elements and subgroups

```@docs
is_conjugate(G::GAPGroup, x::GAPGroupElem, y::GAPGroupElem)
is_conjugate(G::GAPGroup, H::GAPGroup, K::GAPGroup)
representative_action(G::GAPGroup, x::GAPGroupElem, y::GAPGroupElem)
representative_action(G::GAPGroup, H::GAPGroup, K::GAPGroup)
centralizer(G::GAPGroup, x::GAPGroupElem)
centralizer(G::T, H::T) where T <: GAPGroup
normalizer(G::GAPGroup, x::GAPGroupElem)
normalizer(G::T, H::T) where T<:GAPGroup
core(G::T, H::T) where T<:GAPGroup
normal_closure(G::T, H::T) where T<:GAPGroup
```

```@docs
GroupConjClass{T<:GAPGroup, S<:Union{GAPGroupElem,GAPGroup}}
representative(G::GroupConjClass)
acting_group(G::GroupConjClass)
number_conjugacy_classes(G::GAPGroup)
conjugacy_class(G::GAPGroup, g::GAPGroupElem)
conjugacy_class(G::T, g::T) where T<:GAPGroup
conjugacy_classes(G::GAPGroup)
conjugacy_classes_subgroups(G::GAPGroup)
conjugacy_classes_maximal_subgroups(G::GAPGroup)
```


## Cosets (left/right/double)

```@docs
GroupCoset
right_coset(H::GAPGroup, g::GAPGroupElem)
left_coset(H::GAPGroup, g::GAPGroupElem)
is_right(c::GroupCoset)
is_left(c::GroupCoset)
is_bicoset(C::GroupCoset)
acting_domain(C::GroupCoset)
representative(C::GroupCoset)
right_cosets(G::GAPGroup, H::GAPGroup)
left_cosets(G::GAPGroup, H::GAPGroup)
right_transversal(G::T, H::T) where T<: GAPGroup
left_transversal(G::T, H::T) where T<: GAPGroup
GroupDoubleCoset{T <: GAPGroup, S <: GAPGroupElem}
double_coset(G::T, g::GAPGroupElem{T}, H::T) where T<: GAPGroup
double_cosets(G::T, H::T, K::T; check::Bool) where T<: GAPGroup
left_acting_group(C::GroupDoubleCoset)
right_acting_group(C::GroupDoubleCoset)
representative(C::GroupDoubleCoset)
order(C::Union{GroupCoset,GroupDoubleCoset})
Base.rand(C::Union{GroupCoset,GroupDoubleCoset})
intersect(V::AbstractVector{Union{T, GroupCoset, GroupDoubleCoset}}) where T <: GAPGroup
```
