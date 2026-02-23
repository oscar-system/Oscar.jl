```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Nikulin's theory on primitive embeddings

We introduce here the necessary definitions and results which lie behind the
methods about primitive embeddings. Most of the content is taken from
[Nik79](@cite).

## Primitive embeddings

Given an embedding $i\colon S\hookrightarrow T$ of nondegenerate integer
lattices, we call $i$ *primitive* if its cokernel $T/i(S)$ is torsionfree.
There is a left action of the orthogonal group $O(T)$ and a right action
of the orthogonal group $O(S)$ on the set of primitive embeddings as before.
A problem which often arises is to determine representatives for the associated
double cosets.

In [Nik79](@cite), Nikulin gives necessary and sufficient conditions for an
even integer lattice $M$ to embed primitively into an even unimodular lattice
with given invariants (see Theorem 1.12.2 in [Nik79](@cite)). More generally,
the author also provides methods to compute primitive embeddings of any even
lattice into an even lattice of a given genus (see Proposition 1.15.1 in
[Nik79](@cite)). In this proposition, it is explained how to determine
representatives for the double cosets of primitive embeddings as explained in
the first paragraph. It turns out that Nikulin's approach generalizes to the
following.

Let $M$ be an integral integer lattice and let $G$ be a genus of integral
integer lattices (everything nondegenerate). We denote by $L_1, \ldots, L_n$
a complete list of pairwise non-isometric lattices representing all the
isometry classes in $G$. Let $O_M\subset O(M)$ be a subgroup of isometries,
containing the kernel $O^\#(M)$ of the natural map $O(M) \to O(D_M)$,
where $D_M$ denotes the discriminant group of $M$. Then, following Nikulin's
theory on primitive embeddings, there exists an algorithm which returns a
complete set of representatives for the double cosets

```math
O_M\backslash \{M\hookrightarrow L \text{ primitive embedding}\}/O(L)
```
where $L$ runs over all the $L_i$'s.

!!! note
    The lattice $M$ may not admit a primitive embedding in each of the
    $L_i$'s. In particular, if $G$ consists of more than one isometry class,
    Nikulin's approach does not allow to choose into which of the $L_i$'s we
    embed $M$: This would require instead to filter the output using some
    isometry test. However, if $M$ embeds into one of the $L_i$'s, then any
    lattice in the same genus as $M$ also embeds primitively some of the
    $L_j$'s (the nature of Nikulin's approach actually only see local
    invariants of lattices; global invariants such as $O(M)$ are only used
    for a classification purpose).

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
a lattice $L$ in $G$, we construct a lattice $T$ which is unique in its
genus, such that $D_T$ and $D_L(-1)$ are isometric (note that $D_L$ only
depends on $G$) and such that $O(T)\to O(D_T)$ is surjective. We then
classify *primitive extensions* of $M\oplus T$, up to the left action of $O_M$
the left action $O(T)$.

## Primitive extensions

We recall that a *primitive extension* of the orthogonal direct sum of two
integral $\mathbb Z$-lattices $S$ and $T$ is an overlattice $L$ of $M\oplus N$
such that both $S$ and $T$ embed primitively in $L$ (via the natural embeddings
$S,T \to S\oplus T\subseteq L$). Such primitive extensions are obtained, and
classified, by classifying *gluings* between anti-isometric subgroups of the
respective discriminant groups of $S$ and $T$. The corresponding overlattice=
is determined by the graph of such a gluing.

In Proposition 1.5.1 of [Nik79](@cite), Nikulin describes the basis for an
algorithm which answers to the following problem. Given $S$ and $T$ two
integral $\mathbb Z$-lattices and given $O^\#(S)\subset O_S\subset O(S)$ and
$O^\#(T)\subset O_T\subset O(T)$ subgroups of isometries, return a complete
set of representatives for the double cosets in

```math
O_S\backslash \{S\oplus T\subset L \text{ primitive extension}\}/O_T.
```

!!! note
    As in the case of primitive embeddings, since Nikulin's approach focuses on
    working with the discriminant groups of lattices involved, all of our
    classifications are up to the actions of the group of isometries acting
    trivially on said discriminant groups. Therefore, in order to make sense,
    we need to add the action of such groups to our problems.

An algorithm answering to the previous problem, based on Nikulin's work, has
been implemented and it is accessible via the following function.

```@docs
primitive_extensions(::ZZLat, ::ZZLat)
```

## Equivariant primitive extensions

An *equivariant primitive extension* of a pair of $\mathbb Z$-lattices with
isometries $(M, f_M)$ and $(N, f_N)$ is a primitive extension of $M$ and $N$
obtained by gluing two subgroups which are respectively $D_{f_M}$ and
$D_{f_N}$ stable along a glue map which commutes with these two actions.
If such a gluing exists, then the overlattice $L$ of $M\oplus N$ is equipped
with an isometry $f_L$ which preserves both $M$ and $N$, and restricts to $f_M$
and $f_N$ respectively.

Similarly to plain primitive extensions, the work of Nikulin provides all the
ingredients to build an algorithm which computes representatives for all double
cosets of equivariant primitive extensions, for given group actions. Note that
we can actually go a bit further.

If we are given an integral $\mathbb Z$-lattice with isometry $(M, f_M)$ and a
**definite** integral $\mathbb Z$-lattice $N$, we can as well classify primitive
extensions $L$ of $M\oplus N$ equipped with an isometry $f_L$ preserving
$M$ and whose restriction to $M$ is conjugate to $f_M$. This requires an
extra step: Once a representative $L$ of a double coset of primitive extensions
is computed, the algorithm compute representatives for conjugacy classes of
isometries of $N$ (hence requiring $N$ to be definite) which can be extended
by $f_M$ to an isometry of $L$.

The two problems describe above can be answered by calling the following
function.

```@docs
equivariant_primitive_extensions(::Union{ZZLatWithIsom, ZZLat}, ::Union{ZZLatWithIsom, ZZLat})
```

## Admissible equivariant primitive extensions

The following function is a major tool provided by [BH23](@cite). Given
a triple of even $\mathbb Z$-lattices with isometry $((A, a), (B, b), (C, c))$
and two prime numbers $p$ and $q$ (possibly equal), if $(A, B, C)$ is
$p$-admissible, this function returns representatives of isomorphism classes of
equivariant primitive extensions $(A, a)\oplus (B, b)\to (D, d)$ such that the
type of $(D, d^q)$ is equal to the type of $(C, c)$
(see [`type(::ZZLatWithIsom)`](@ref)).

```@docs
admissible_equivariant_primitive_extensions(::ZZLatWithIsom, ::ZZLatWithIsom, ::ZZLatWithIsom, ::Int, ::Int)
```
