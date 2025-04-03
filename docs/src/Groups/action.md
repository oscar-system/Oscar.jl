```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Group actions

A *group action* of a group $G$ on a set $\Omega$ (from the right) is defined by
a map $\mu:\Omega\times G\to \Omega$ that satisfies the compatibility conditions
$\mu(\mu(x,g),h) = \mu(x, gh)$ and $\mu(x, 1_G) = x$ for all $x\in\Omega$.

The maps $\mu$ are implemented as functions that take two arguments, an element
$x$ of $\Omega$ and a group element $g$, and return the image of $x$ under $g$.

In many cases, a natural action is given by the types of the elements in $\Omega$
and in $G$.
For example permutation groups act on positive integers by just applying
the permutations.
In such situations, the function `^` can be used as action function,
and `^` is taken as the default whenever no other function is prescribed.

However, the action is not always determined by the types of the involved
objects.
For example, permutations can act on vectors of positive integers by
applying the permutations pointwise, or by permuting the entries;
matrices can act on vectors by multiplying the vector with the matrix,
or by multiplying the inverse of the matrix with the vector;
and of course one can construct new custom actions in situations where
default actions are already available.

Thus it is in general necessary to specify the action function explicitly,
see the following sections.


## Common actions of group elements

```@docs
on_tuples
on_sets
permuted
on_indeterminates
on_lines
on_echelon_form_mats
on_subgroups
```


## G-Sets

The idea behind G-sets is to have objects that encode the permutation action
induced by a group (that need not be a permutation group) on a given set.
A G-set provides an explicit bijection between the elements of the set and
the corresponding set of positive integers on which the induced permutation
group acts,
see [`action_homomorphism(Omega::GSetByElements{T}) where T<:GAPGroup`](@ref).
Note that the explicit elements of a G-set `Omega` can be obtained using
`collect(Omega)`.

```@docs
gset(G::Union{GAPGroup, FinGenAbGroup}, fun::Function, Omega)
permutation
acting_group(Omega::GSetByElements)
action_function(Omega::GSetByElements)
action_homomorphism(Omega::GSetByElements{T}) where T<:GAPGroup
is_conjugate(Omega::GSet, omega1, omega2)
is_conjugate_with_data(Omega::GSet, omega1, omega2)
orbit(Omega::GSetByElements{<:GAPGroup, S}, omega::S) where S
orbit(G::PermGroup, omega)
orbits(Omega::T) where T <: GSetByElements{TG} where TG <: GAPGroup
is_transitive(Omega::GSet)
transitivity(Omega::GSet)
rank_action(Omega::GSet)
is_primitive(Omega::GSet)
is_regular(Omega::GSet)
is_semiregular(Omega::GSet)
```

## Block systems of a G-set

If we have a G-set $\Omega$, a *block system* of $\Omega$ is a partition that is invariant under the action of the associated group.
The group action on $\Omega$ induces a natural action on such a partition.

When calling these methods with a `GSet` as the argument, we require that the group action is transitive.
The blocks are returned as Julia `Set` objects.
Note that this is in contrast to the return type when calling the methods with a `PermGroup` as the argument, in which case the blocks are sorted vectors of integers.

```@docs
blocks(Omega::GSet)
maximal_blocks(Omega::GSet)
minimal_block_reps(Omega::GSet)
all_blocks(Omega::GSet)
```


## Stabilizers

```@docs
stabilizer(G::GAPGroup, pnt::Any, actfun::Function)
stabilizer(Omega::GSet)
```
