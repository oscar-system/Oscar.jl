```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Subgroups](@id subgroups)

The following functions are available in OSCAR for subgroup properties:

```@docs
sub(G::GAPGroup, gens::AbstractVector{<:GAPGroupElem}; check::Bool = true)
is_subset(H::GAPGroup, G::GAPGroup)
is_subgroup(H::GAPGroup, G::GAPGroup)
embedding(H::GAPGroup, G::GAPGroup)
index(G::GAPGroup, H::GAPGroup)
is_maximal_subgroup(H::GAPGroup, G::GAPGroup; check::Bool = true)
is_normalized_by(H::GAPGroup , G::GAPGroup)
is_normal_subgroup(H::GAPGroup, G::GAPGroup)
is_characteristic_subgroup(H::GAPGroup, G::GAPGroup; check::Bool = true)
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
socle
solvable_radical
pcore(G::GAPGroup, p::IntegerUnion)
intersect(::T, V::T...) where T<:GAPGroup
```

The following functions return a vector of subgroups.

```@docs
normal_subgroups
maximal_normal_subgroups
minimal_normal_subgroups
characteristic_subgroups
derived_series(G::GAPGroup)
sylow_system
hall_system
complement_system
chief_series
composition_series
jennings_series
p_central_series
lower_central_series(G::GAPGroup)
upper_central_series
```

!!! note
    When a function returns a vector of subgroups,
    the output consists in the subgroups only;
    the embeddings are not returned as well.
    To get the embedding homomorphism of the subgroup `H` in `G`,
    one can type `embedding(H, G)`.


The following functions return an iterator of subgroups.
Usually it is more efficient to work with (representatives of) the
underlying conjugacy classes of subgroups instead.

```@docs
complements(G::GAPGroup, N::GAPGroup)
hall_subgroups
low_index_subgroups
maximal_subgroups
subgroups(G::GAPGroup)
```

## Conjugation action of elements and subgroups

```@docs
is_conjugate(G::GAPGroup, x::GAPGroupElem, y::GAPGroupElem)
is_conjugate(G::GAPGroup, H::GAPGroup, K::GAPGroup)
is_conjugate_with_data(G::GAPGroup, x::GAPGroupElem, y::GAPGroupElem)
is_conjugate_with_data(G::GAPGroup, H::GAPGroup, K::GAPGroup)
centralizer(G::GAPGroup, x::GAPGroupElem)
centralizer(G::GAPGroup, H::GAPGroup)
normalizer(G::GAPGroup, x::GAPGroupElem)
normalizer(G::GAPGroup, H::GAPGroup)
core(G::GAPGroup, H::GAPGroup)
normal_closure(G::GAPGroup, H::GAPGroup)
```

```@docs
GroupConjClass{T<:GAPGroup, S<:Union{GAPGroupElem,GAPGroup}}
representative(G::GroupConjClass)
acting_group(G::GroupConjClass)
number_of_conjugacy_classes(G::GAPGroup)
conjugacy_class(G::GAPGroup, g::GAPGroupElem)
conjugacy_class(G::GAPGroup, H::GAPGroup)
conjugacy_classes(G::GAPGroup)
complement_classes
hall_subgroup_classes
low_index_subgroup_classes
maximal_subgroup_classes(G::GAPGroup)
subgroup_classes(G::GAPGroup)
```


## Cosets (left/right/double)

```@docs
GroupCoset
group(C::GroupCoset)
acting_group(C::GroupCoset)
representative(C::GroupCoset)
right_coset(H::GAPGroup, g::GAPGroupElem)
left_coset(H::GAPGroup, g::GAPGroupElem)
is_right(C::GroupCoset)
is_left(C::GroupCoset)
right_cosets(G::GAPGroup, H::GAPGroup; check::Bool=true)
left_cosets(G::GAPGroup, H::GAPGroup; check::Bool=true)
right_transversal(G::T1, H::T2; check::Bool=true) where T1 <: GAPGroup where T2 <: GAPGroup
left_transversal(G::T1, H::T2; check::Bool=true) where T1 <: GAPGroup where T2 <: GAPGroup
is_bicoset(C::GroupCoset)
```

```@docs
GroupDoubleCoset{T <: GAPGroup, S <: GAPGroupElem}
group(C::GroupDoubleCoset)
left_acting_group(C::GroupDoubleCoset)
right_acting_group(C::GroupDoubleCoset)
representative(C::GroupDoubleCoset)
double_coset(G::GAPGroup, g::GAPGroupElem, H::GAPGroup)
double_cosets(G::T, H::GAPGroup, K::GAPGroup; check::Bool=true) where T <: GAPGroup
```

```@docs
order(C::Union{GroupCoset,GroupDoubleCoset})
Base.rand(C::Union{GroupCoset,GroupDoubleCoset})
```
