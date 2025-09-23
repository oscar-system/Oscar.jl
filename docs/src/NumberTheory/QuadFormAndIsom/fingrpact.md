```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Collections of isometries

In this section, we review methods related to collection of isometries for
rational quadratic spaces and integer lattices.

## Change of basis representation

Given an integer lattice $L$ and a collection $F$ of isometries of $L$,
there are two canonical ways of representing the isometries in $F$ by mean
of matrices. The first way is to represent the isometries in $F$ in the
standard basis of the ambient quadratic space of $L$. The second way is to
represent such isometries in the given fixed basis of $L$.

Given one representation of $F$, by mean of matrices, the following two methods
allow to switch to another one.

```@docs
extend_to_ambient_space
representation_in_lattice_coordinates
```

## Invariant and coinvariant lattices

We provide the possibility to compute invariant and coinvariant
sublattices given an orthogonal representation `G` in matrix form of a
collection of isometries of a given lattice `L`:

```@docs
coinvariant_lattice(::ZZLat, ::MatrixGroup)
invariant_lattice(::ZZLat, ::MatrixGroup)
invariant_coinvariant_pair(::ZZLat, ::Union{QQMatrix, Vector{QQMatrix}, MatrixGroup})
```

## Special isometries

Given an integer lattice $L$, any isometry $f$ of $L$ has determinant equal to
$\pm 1$. We define the *special_orthogonal_group* $SO(L)$ of $L$ to be the
normal subgroup of $O(L)$ consisting of isometries of determinant $+1$.

```@docs
special_orthogonal_group(::ZZLat)
```

Similarly, given any finite group $G$ of isometries of $L$, one can compute
the special subgroup $G\cap SO(L)$ of $G$.

```@docs
special_subgroup(::ZZLat, ::MatrixGroup)
```

## Stable isometries

Given an integral integer lattice $L$, any isometry $f$ of $L$ acts on the
discriminant group $D_L$ of $L$. We define the *stable_orthogonal_group*
$O^\#(L)$ of $L$ to be the normal subgroup of $O(L)$ consisting of isometries
acting trivially on $D_L$.

```@docs
stable_orthogonal_group(::ZZLat)
```

Similarly, given any finite group $G$ of isometries of $L$, one can compute
the stable subgroup $G\cap O^\#(L)$ of $G$.

```@docs
stable_subgroup(::ZZLat, ::MatrixGroup)
```

## Stabilizers

```@docs
stabilizer_discriminant_subgroup
stabilizer_in_diagonal_action
maximal_extension(::ZZLat, ::ZZLat, ::MatrixGroup)
stabilizer_in_orthogonal_group
pointwise_stabilizer_in_orthogonal_group
setwise_stabilizer_in_orthogonal_group
pointwise_stabilizer_orthogonal_complement_in_orthogonal_group
```

## Saturation

Given an integer lattice $L$ and two groups of isometries $H \leq G$ of $L$,
we define the *saturation* of $H$ in $G$ to be the kernel of group homomorphism

```math
G(L^H) \to O(L^H)
```
where $G(L^H)$ denotes the stabilizer of the invariant sublattice $L^H$ in $G$.
In other words, the saturation of $H$ in $G$ is the pointwise stabilizer of
$L^H$ in $G$. In the case where the group $G$ is finite, one can actually
explictly compute such a group:

```@docs
saturation(::ZZLat, ::MatrixGroup, ::MatrixGroup)
```
In general, if the coinvariant sublattice $L_H$ of the group $H$ is definite or
of rank 2, one can compute the saturation of $H$ inside the orthogonal group
$O(L)$ of $L$.

```@docs
saturation(::ZZLat, ::MatrixGroup)
```
Finally, whenever the coinvariant lattice $L_H$ of $H$ is definite, it is
possible to decide whether $H$ is saturated in the given groups as above.

```@docs
is_saturated_with_saturation
```

## Isometry checks

Given a rational quadratic space or an integer lattice, one can ask whether
a collection of matrices with rational entries consists of isometries of the
given object. The following six functions are not exported, they can be used
for input checks in one's functions.

```@docs
is_isometry(::Hecke.QuadSpace, ::QQMatrix)
is_isometry(::ZZLat, ::QQMatrix)
is_isometry_list(::Hecke.QuadSpace, ::Vector{QQMatrix})
is_isometry_list(::ZZLat, ::Vector{QQMatrix})
is_isometry_group(::Hecke.QuadSpace, ::MatrixGroup)
is_isometry_group(::ZZLat, ::MatrixGroup)
```

Given an integer lattice with isometry $(L, f)$, one can also ask whether
the isometry $f$ satisfies some properties.

```@docs
is_stable_isometry(::ZZLatWithIsom)
is_special_isometry(::ZZLatWithIsom)
```

