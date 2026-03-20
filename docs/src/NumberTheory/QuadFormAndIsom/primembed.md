```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Nikulin's theory on primitive embeddings

We introduce here the necessary definitions and results which lie behind the
methods about primitive embeddings. Most of the content is taken from the first
section of [Nik79](@cite). This chapter only deals with nondegenerate
$\mathbb{Z}$-lattices, which we simply call lattices.

## Primitive embeddings

Given an embedding $i\colon M\hookrightarrow L$ of lattices, we call
$i$ **primitive** if its cokernel $L/i(M)$ is torsionfree.
There is a left action of the orthogonal group $O(L)$ and a right action
of the orthogonal group $O(M)$ on the set of primitive embeddings as before.
A problem which often arises is to determine representatives for the associated
double cosets.

In [Nik79](@cite), Nikulin gives necessary and sufficient conditions for an
even lattice $M$ to embed primitively into an even unimodular lattice with
given invariants (see Theorem 1.12.2 in [Nik79](@cite)). More generally,
the author also provides methods to compute primitive embeddings
of any even lattice into an even lattice of a given genus (see Proposition
1.15.1 in [Nik79](@cite)). In this proposition, it is explained how to
determine representatives for the double cosets of primitive embeddings as
explained in the first paragraph. It turns out that Nikulin's approach
generalizes to the following.

Let $M$ be an integral lattice and let $G$ be a genus of integral
lattices. We denote by $L_1,\ldots, L_n$ a complete list of
pairwise non-isometric lattices representing all the isometry classes in $G$.
Let $O_M\subset O(M)$ be a subgroup of isometries, containing the kernel
$O^\#(M)$ of the natural map $O(M) \to O(D_M)$, where $D_M$ denotes the
discriminant group of $M$. Then, following Nikulin's theory on primitive
embeddings, there exists an algorithm which returns a complete set of
representatives for the double cosets

```math
O(L)\backslash \{M\hookrightarrow L \text{ primitive embedding}\}/O_M
```
where $L$ runs over all the $L_i$'s.

!!! note
    The lattice $M$ may not admit a primitive embedding in each of the
    $L_i$'s. In particular, if $G$ consists of more than one isometry class,
    Nikulin's approach does not allow to choose into which of the $L_i$'s the
    lattice $M$ embeds: This would require instead to filter the output using
    some isometry test. However, if $M$ embeds into one of the $L_i$'s, then
    any lattice in the same genus as $M$ also embeds primitively into some of
    the $L_j$'s (the nature of Nikulin's approach is local; global invariants
    such as $O(M)$ are only used for the classification purpose at hand).

Such an algorithm has been designed and implemented as part of this package,
and it can be accessed via the following function.

```@docs
primitive_embeddings(::ZZGenus, ::ZZLat)
```

!!! warning
    The algorithm tends to be slow for large rank or determinant. This has to
    do with computations of orbits and stabilizers of subgroups in finite
    abelian groups.

The idea behind the algorithm is a generalization of the proof of
Proposition 1.15.1 of [Nik79](@cite). In order to embed $M$ primitively in
a lattice $L$ in $G$, we construct a lattice $N$ which is unique in its
genus, such that $D_N$ and $D_L(-1)$ are isometric (note that $D_L$ only
depends on $G$) and such that $O(N)\to O(D_N)$ is surjective. We then
classify *primitive extensions* of $M$ and $N$, up to the right action of
$O_M$ and the left action of $O(N)$.

## Primitive extensions

A **primitive extension** of a finite collection of lattices
$\{M_i\}_{i=1,\ldots,n}$ is an overlattice $L$ of $\bigoplus_{i=1}^nM_i$
such that each $M_i$ is primitive in $L$.
Given two integral lattices $M$ and $N$, their primitive extensions are in
bijection with anti-isometries between subgroups of their respective
discriminant groups. The corresponding overlattice is determined by the graph
of such an anti-isometry, also known as a **glue map**.

In Proposition 1.5.1 of [Nik79](@cite), Nikulin describes the basis for an
algorithm which answers to the following problem. Given $M$ and $N$ two
integral lattices and given $O^\#(M)\subset O_M\subset O(M)$ and
$O^\#(N)\subset O_N\subset O(N)$ subgroups of isometries, return a complete
set of representatives for the double cosets in

```math
O_N\backslash \{M\oplus N\subset L \text{ primitive extension}\}/O_M.
```

!!! note
    As in the case of primitive embeddings, since Nikulin's approach focuses on
    working with the discriminant groups of lattices involved, all of the
    classifications are up to the actions of the group of isometries acting
    trivially on said discriminant groups. Therefore, in order to make sense,
    we need to add the action of such groups to the problems.

An algorithm answering to the previous problem, based on Nikulin's work, is
available via the following function.

```@docs
primitive_extensions(::ZZLat, ::ZZLat)
```

## Equivariant primitive extensions

An **equivariant primitive extension** of a finite family of pairs of
lattices with isometry $\{(M_i, f_i)\}_{i=1,\ldots,n}$ is a primitive extension
$L$ of $\{M_i\}_{i=1,\ldots,n}$ preserved by the diagonal action of the $f_i$'s.
Given two pairs $(M, f_M)$ and $(N, f_N)$ where $M$ and $N$ are integral,
their equivariant primitive extensions are in bijection with glue maps $\gamma$
between subgroups $H_M\leq D_M$ and $H_N\leq D_N$ which are stable under the
respective actions of $f_M$ and $f_N$, and for which the following holds:
```math
\gamma\circ D_{f_M}|_{H_M} = D_{f_N}|_{H_N}\circ \gamma.
```
In fact, given any such glue map, the corresponding overlattice $L$ of
$M\oplus N$ is equipped with an isometry $f_L$ which preserves both $M$ and
$N$, and restricts to $f_M$ and $f_N$ respectively.

Similarly to plain primitive extensions, the work of Nikulin provides all the
ingredients to build an algorithm which computes representatives for all double
cosets of equivariant primitive extensions, for given group actions. Note that
we can actually go a bit further.

If we are given an integral lattice with isometry $(M, f_M)$ and a *definite*
integral lattice $N$, one can as well classify primitive extensions $L$ of $M$
and $N$ which are equipped with an isometry $f_L$ preserving $M$ and whose
restriction to $M$ is conjugate to $f_M$. This requires an extra step: Once a
representative $L$ of a double coset of primitive extensions is computed, the
algorithm computes representatives for conjugacy classes of isometries of $N$
(hence requiring $N$ to be definite) which can be extended by $f_M$ to an
isometry of $L$.

The two problems described above can be solved by calling the following
function.

```@docs
equivariant_primitive_extensions(::Union{ZZLatWithIsom, ZZLat}, ::Union{ZZLatWithIsom, ZZLat})
```

## Admissible equivariant primitive extensions

The following function is a major tool provided by [BH23](@cite). Given
a triple of even lattices with isometry $((A, a), (B, b), (C, c))$ and two
prime numbers $p$ and $q$ (possibly equal), if $(A, B, C)$ is $p$-admissible,
this function returns representatives of isomorphism classes of equivariant
primitive extensions $(A, a)\oplus (B, b)\to (D, d)$ such that the type of
$(D, d^q)$ is equal to the type of $(C, c)$
(see [`type(::ZZLatWithIsom)`](@ref)).

```@docs
admissible_equivariant_primitive_extensions(::ZZLatWithIsom, ::ZZLatWithIsom, ::ZZLatWithIsom, ::Int, ::Int)
```
