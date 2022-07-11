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
Pages = ["permgroup.md"]
```

# Permutation groups

Permutation groups can be defined as symmetric groups, alternating groups or their subgroups.

```@docs
PermGroup
PermGroupElem
symmetric_group
is_natural_symmetric_group(G::GAPGroup)
is_isomorphic_with_symmetric_group(G::GAPGroup)
alternating_group
is_natural_alternating_group(G::GAPGroup)
is_isomorphic_with_alternating_group(G::GAPGroup)
```

In Oscar, every permutation group has a degree `n`, that corresponds to the size of the set on which `G` acts.

```@docs
degree(x::PermGroup)
```

## Permutations

Permutations in Oscar are displayed as products of disjoint cycles, as in GAP. An explicit permutation can be built using the functions `perm`, `gap_perm`, `cperm`, or `@perm`.

```@docs
perm
gap_perm
cperm
@perm
```

Every permutation has always a permutation group as a parent. Two permutations coincide if, and only if, they move the same points and their parent groups have the same degree.
```jldoctest
julia> G=symmetric_group(5);

julia> A=alternating_group(5);

julia> x=cperm(G,[1,2,3]);

julia> y=cperm(A,[1,2,3]);

julia> z=cperm([1,2,3]); parent(z)
Sym( [ 1 .. 3 ] )

julia> x==y
true

julia> x==z
false
```
In the example above, `x` and `y` are equal because both act on a set of cardinality `5`, while `x` and `z` are different because `x` belongs to `Sym(5)` and `z` belongs to `Sym(3)`.

One can also create permutation group elements with the following syntax
```jldoctest
julia> G=symmetric_group(6)
Sym( [ 1 .. 6 ] )

julia> G([2,4,6,1,3,5])
(1,2,4)(3,6,5)
```

Here `G` is the parent group. 
```jldoctest
julia> G=symmetric_group(6)
Sym( [ 1 .. 6 ] )

julia> x=G([2,4,6,1,3,5])
(1,2,4)(3,6,5)

julia> y=perm(G,[2,4,6,1,3,5])
(1,2,4)(3,6,5)

julia> x==y
true
```
If `G` is a group and `x` is a permutation,
`G(x)` returns a permutation `x` with parent `G`;
an exception is thrown if `x` does not embed into `G`.
```@repl oscar
G=symmetric_group(5);
x=cperm([1,2,3]);
y=G(x);
parent(x)
parent(y)
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

## Operations for permutation groups

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
