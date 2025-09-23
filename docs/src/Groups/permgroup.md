```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Permutation groups

## Constructing permutation groups

Permutation groups can be defined as symmetric groups, alternating groups or their subgroups.

```@docs
PermGroup
PermGroupElem
symmetric_group
alternating_group
permutation_group
@permutation_group
```

```@docs
projective_general_linear_group
projective_special_linear_group
projective_symplectic_group
projective_orthogonal_group
projective_special_orthogonal_group
projective_omega_group
projective_unitary_group
projective_special_unitary_group
```

## Operations for permutation groups

All operations, properties and attributes for general groups
described in the previous sections are supported for permutation
groups. In addition there are some specific to permutation groups.

In OSCAR, every permutation group has a degree `n`, that corresponds to the size of the set on which `G` acts.

```@docs
degree(x::PermGroup)
smaller_degree_permutation_representation(G::PermGroup)
```


The following functions deal with the natural action of a given permutation group $G$.
```@docs
is_transitive(G::PermGroup, L::AbstractVector{Int} = 1:degree(G))
transitivity(G::PermGroup, L::AbstractVector{Int} = 1:degree(G))
is_primitive(G::PermGroup, L::AbstractVector{Int} = 1:degree(G))
is_regular(G::PermGroup, L::AbstractVector{Int} = 1:degree(G))
is_semiregular(G::PermGroup, L::AbstractVector{Int} = 1:degree(G))
rank_action(G::PermGroup, L::AbstractVector{Int} = 1:degree(G))
blocks(G::PermGroup, L::AbstractVector{Int} = moved_points(G))
maximal_blocks(G::PermGroup, L::AbstractVector{Int} = moved_points(G))
minimal_block_reps(G::PermGroup, L::AbstractVector{Int} = moved_points(G))
all_blocks(G::PermGroup)
```

The following functions allow efficiently "recognizing" certain permutation groups as alternating or symmetric groups.

```@docs
is_natural_symmetric_group(G::GAPGroup)
is_isomorphic_to_symmetric_group(G::GAPGroup)
is_natural_alternating_group(G::GAPGroup)
is_isomorphic_to_alternating_group(G::GAPGroup)
```

## Permutations

Permutations in OSCAR are displayed as products of disjoint cycles, as in GAP. An explicit permutation can be built using the functions `perm`, `cperm`, or `@perm`.

```@docs
perm
cperm
@perm
```

The function `Vector{T}` works in the opposite way with respect to `perm`:
```@docs
Vector(x::PermGroupElem, n::Int = x.parent.deg)
```

## Operations on permutations

```@docs
sign(g::PermGroupElem)
isodd(g::PermGroupElem)
iseven(g::PermGroupElem)
cycle_structure(g::PermGroupElem)
cycles(g::PermGroupElem)
```

## Permutations as functions
A permutation can be viewed as a function on the set `{1,...,n}`, hence it can be evaluated on integers.

!!! note
    The multiplication between permutations works from the left to the right. So, if `x` and `y` are permutations and `n` is an integer, then `(x*y)(n) = (y(x(n))`, NOT `x(y(n))`.
    This works also if the argument is not in the range `1:n`; in such a case, the output coincides with the input.

```jldoctest
julia> x = cperm([1,2,3,4,5]);

julia> x(2)
3

julia> x(6)
6
```

## Moved points of permutations and permutation groups
```@docs
smallest_moved_point
largest_moved_point
moved_points
number_of_moved_points
fixed_points
number_of_fixed_points
```

## Cycle structures

For a permutation, its cycle structure [`cycle_structure`](@ref)
determines the degree, order, number of moved points, sign.

```@docs
degree(::CycleType)
iseven(::CycleType)
isodd(::CycleType)
order(::Type{T}, c::CycleType) where T
sign(::CycleType)
```

```@docs
cycle_structures(G::PermGroup)
```
