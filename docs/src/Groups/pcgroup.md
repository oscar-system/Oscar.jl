# Polycyclic groups

## Introduction

A *polycyclic group* is a group $G$ which has a subnormal series
$G = G_1 > G_2 > \cdots > G_n > G_{n+1} = \{ 1 \}$
such that the factor groups $G_i/G_{i+1}$ are cyclic.

Polycyclic groups are solvable, and every finite solvable group is polycyclic.

Each polycyclic group has a finite presentation with the following structure:
Choose elements $g_i \in G_i \setminus G_{i+1}$, for $1 \leq i \leq n$,
such that $G_i = \langle G_{i+1}, g_i \rangle$ holds,
and let $r_i$ be the index of $G_{i+1}$ in $G_i$ (a positive integer or
infinity), the *relative order* of $g_i$ w.r.t. the given subnormal series.
Then each element of $G$ has a unique expression of the form
$g_1^{e_1} g_2^{e_2} \cdots g_n^{e_n}$ where the exponent $e_i$ is in
the set $\{ 0, 1, \ldots, r_i-1 \}$ if $r_i$ is finite, and $e_i$ is an
integer otherwise; this expression is called the *normal form*
of the group element.

The group $G$ has a presentation with generators $g_1, g_2, \ldots, g_n$
and defining relations of the following form,
where we define $I$ as the set of indices $i$ for which $r_i$ is a positive integer.

```math
   \begin{array}{lcll}
      g_i^{r_i} & = & g_{i+1}^{a(i,i,i+1)} \cdots g_n^{a(i,i,n)}, &
         1 \leq i \leq n, i \in I \\
      g_i^{-1} g_j g_i & = & g_{i+1}^{a(i,j,i+1)} \cdots g_n^{a(i,j,n)}, &
         1 \leq i < j \leq n.
   \end{array}
```

For infinite groups, we need additionally

```math
   \begin{array}{lcll}
      g_i^{-1} g_j^{-1} g_i & = & g_{i+1}^{b(i,j,i+1)} \cdots g_n^{b(i,j,n)}, &
         1 \leq i < j \leq n, j \not\in I \\
      g_i g_j g_i^{-1} & = & g_{i+1}^{c(i,j,i+1)} \cdots g_n^{c(i,j,n)}, &
         1 \leq i < j \leq n, i \not\in I \\
      g_i g_j^{-1} g_i^{-1} & = & g_{i+1}^{d(i,j,i+1)} \cdots g_n^{d(i,j,n)}, &
         1 \leq i < j \leq n, i, j, \not\in I.
   \end{array}
```

Here we assume that the right hand sides are words in normal form.

This presentation describes a confluent rewriting system,
thus it admits the computation of the normal form for each group element.
Inverses and products of elements given by their normal forms
can be efficiently written in normal form.

Groups which are given by a polycyclic presentation are called
*polycyclicly presented* groups.
The rest of this section is about these groups.

Polycyclicly presented groups in Oscar have the type [`PcGroup`](@ref),
their elements have the type [`PcGroupElem`](@ref).
Analogous to the situation with finitely presented groups and their subgroups,
there are the types [`SubPcGroup`](@ref) for subgroups
and [`SubPcGroupElem`](@ref) for their elements.

## Basic Creation

One can write down a polycyclic presentation by hand, using a so-called
[collector](@ref "Collectors for polycyclicly presented groups"),
and then creating a group with this presentation.

```jldoctest
julia> c = collector(2, Int);

julia> set_relative_orders!(c, [2, 3])

julia> set_conjugate!(c, 2, 1, [2 => 2])

julia> gg = pc_group(c)
Pc group of order 6

julia> describe(gg)
"S3"
```

Alternatively, one can take a polycyclic group,
and let Oscar compute a pc presentation for it.

```jldoctest
julia> g = symmetric_group(4)
Sym(4)

julia> iso = isomorphism(PcGroup, g);

julia> h = codomain(iso)
Pc group of order 24
```

For certain series of groups, one can
[create the members as pc groups](@ref "Series of polycyclic groups").

```jldoctest
julia> dihedral_group(8)
Pc group of order 8
```

And the groups from [the small groups library](@ref "Groups of small order")
are represented by pc groups whenever they are solvable.

```jldoctest
julia> small_group(24, 12)
Pc group of order 24
```

## Functions for elements of (subgroups of) pc groups

```@docs
letters(g::Union{PcGroupElem, SubPcGroupElem})
syllables(g::Union{PcGroupElem, SubPcGroupElem})
map_word(g::Union{PcGroupElem, SubPcGroupElem}, genimgs::Vector; genimgs_inv::Vector = Vector(undef, length(genimgs)), init = nothing)
```

## Functions for (subgroups of) pc groups

```@docs
relators(G::PcGroup)
```

The function
[`full_group(G::T) where T <: Union{SubFPGroup, SubPcGroup}`](@ref)
for (subgroups of) finitely presented groups is applicable to
(subgroups of) polycyclicly presented groups as well.

## Series of polycyclic groups

The following functions can be used to create polycyclicly presented groups
from certain series of groups.

(In fact, one can request also other types for the results,
but `PcGroup` is the default type in all cases except `abelian_group`.)

```@docs
abelian_group(::Type{T}, v::Vector{Int}) where T <: GAPGroup
elementary_abelian_group
extraspecial_group
cyclic_group
dihedral_group
quaternion_group
```

## Collectors for polycyclicly presented groups

The following functions can be used to enter polycyclic presentations
by hand, or to create new such presentations from given ones.

```@docs
collector(n::Int, ::Type{T} = ZZRingElem) where T <: IntegerUnion
set_relative_order!(c::Collector{T}, i::Int, relord::T) where T <: IntegerUnion
get_relative_order(c::Collector{T}, i::Int) where T <: IntegerUnion
set_relative_orders!(c::Collector{T}, relords::Vector{T}) where T <: IntegerUnion
get_relative_orders(c::Collector{T}) where T <: IntegerUnion
set_power!(c::Collector{T}, i::Int, rhs::Vector{Pair{Int, T}}) where T <: IntegerUnion
get_power(c::Collector{T}, i::Int) where T <: IntegerUnion
set_conjugate!(c::Collector{T}, j::Int, i::Int, rhs::Vector{Pair{Int, T}}) where T <: IntegerUnion
get_conjugate(c::Collector{T}, j::Int, i::Int) where T <: IntegerUnion
set_commutator!(c::Collector{T}, j::Int, i::Int, rhs::Vector{Pair{Int, T}}) where T <: IntegerUnion
pc_group(c::GAP_Collector)
collector(::Type{T}, G::PcGroup) where T <: IntegerUnion
```

## Technicalities

```@docs
PcGroup
PcGroupElem
SubPcGroup
SubPcGroupElem
```
