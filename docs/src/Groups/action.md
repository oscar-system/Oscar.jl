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
orbit(Omega::GSetByElements{<:GAPGroup}, omega::T) where T
orbit(G::PermGroup, omega)
orbits(Omega::T) where T <: GSetByElements{TG} where TG <: GAPGroup
```


## Stabilizers

```@docs
stabilizer(G::Oscar.GAPGroup, pnt::Any, actfun::Function)
```
