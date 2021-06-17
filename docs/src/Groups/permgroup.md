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
alternating_group
```

In Oscar, every permutation group has a degree `n`, that corresponds to the size of the set on which `G` acts.

```@docs
degree(x::PermGroup)
```

!!! note
    The degree of a group of permutations is not necessarily equal to the largest moved point of the group `G`. For example, the trivial subgroup of `symmetric_group(n)` has degree `n` even though it fixes `n`.

## Permutations

Permutations in Oscar are displayed as products of disjoint cycles, as in GAP. An explicit permutation can be built using the functions `perm`, `gap_perm` and `cperm`.

```@docs
perm
gap_perm
cperm
```

  **Example:**
```jldoctest
julia> perm(symmetric_group(6),[2,4,6,1,3,5])
(1,2,4)(3,6,5)

julia> cperm([1,2,3],4:7)
(1,2,3)(4,5,6,7)

julia> cperm([1,2],[2,3])
(1,3,2)
```

At the moment, the input vectors of the function `cperm` need not to be disjoint.

!!! warning
    If the function `perm` is evaluated in a vector of integers without specifying the group `G`, then the returned value is an element of the AbstractAlgebra.jl type `Perm{Int}`. For this reason, if one wants a permutation of type `GAPGroupElem{PermGroup}` without specifying a parent, one has to use the function `gap_perm`.


Every permutation has always a permutation group as a parent. Two permutations coincide if, and only if, they move the same points and their parent groups have the same degree.
```jldoctest
julia> G=symmetric_group(5);

julia> A=alternating_group(5);

julia> x=cperm(G,[1,2,3]);

julia> y=cperm(A,[1,2,3]);

julia> z=cperm([1,2,3]);

julia> x==y
true

julia> x==z
false
```
In the example above, `x` and `y` are equal because both act on a set of cardinality `5`, while `x` and `z` are different because `x` belongs to `Sym(5)` and `z` belongs to `Sym(2)`.


If `G` is a group and `x` is a permutation, it is possible to set `G` as parent of `x` simply typing `G(x)`. This returns the permutation `x` as element of `G` (or ERROR if `x` does not embed into `G`).
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


## Permutations as functions
A permutation can be viewed as a function on the set `{1,...,n}`, hence it can be evaluated on integers.

!!! note
    The multiplication between permutations works from the left to the right. So, if `x` and `y` are permutations and `n` is an integer, then `(x*y)(n) = (y(x(n))`, NOT `x(y(n))`.

```jldoctest
julia> x = cperm([1,2,3,4,5]);

julia> x(2)
3
```
This works also if the argument is not in the range `1:n`; in such a case, the output coincides with the input.
