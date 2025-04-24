```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

We introduce here the necessary definitions and results which lie behind the
methods about primitive embeddings. Most of the content is taken from
[Nik79](@cite).

# Nikulin's theory on primitive embeddings

## Primitive embeddings

Given an embedding $i\colon S\hookrightarrow T$ of non-degenerate integral
integer lattices, we call $i$ *primitive* if its cokernel $T/i(S)$ is torsion
free. Two primitive embeddings $i_1\colon S\hookrightarrow M_1$ and
$i_2\colon S \hookrightarrow M_2$ of $S$ into two lattices $M_1$ and $M_2$ are
called *isomorphic* if there exists an isometry $M_1 \to M_2$ which restricts to
the identity of $S$. Moreover, if there exists an isometry between $M_1$ and
$M_2$ which maps $S$ to itself (not necessarily identically), we say that $i_1$
and $i_2$ defines *isomorphic primitive sublattices* [Nik79](@cite).

In his paper, V. V. Nikulin gives necessary and sufficient condition for an even
integral lattice $M$ to embed primitively into an even unimodular lattice with
given invariants (see Theorem 1.12.2 in [Nik79](@cite)). More generally, the
author also provides methods to compute primitive embeddings of any even lattice
into an even lattice of a given genus (see Proposition 1.15.1 in [Nik79](@cite)).
In the latter proposition, it is explained how to classify such embeddings as
isomorphic embeddings or as isomorphic sublattices. Moreover, with enough care,
one can generalize the previous results for embeddings in odd lattices.

A general method to compute primitive embeddings between integral lattices
can be algorithmically implemented, however it tends to be slow and inefficient
in general for large rank or determinant. But, in the case where the
discriminant groups are (elementary) $p$-groups, the method can be made
more efficient.

We provide 4 kinds of output:
* A boolean, which only returns whether there exists a primitive embedding;
* A single primitive embedding as soon as the algorithm computes one;
* A list of representatives of isomorphism classes of primitive embeddings;
* A list of representatives of isomorphism classes of primitive sublattices.

```@docs
primitive_embeddings(::ZZLat, ::ZZLat)
```

Note that the previous two functions require the first lattice of the input to be
unique in its genus. Otherwise, one can specify a genus, or its invariants, as a
first input:

```@docs
primitive_embeddings(::ZZGenus, ::ZZLat)
primitive_embeddings(::TorQuadModule, ::Tuple{Int, Int}, ::ZZLat)
```

In order to compute such primitive embeddings of a lattice $M$ into a lattice
$L$, we follow the proof of Proposition 1.15.1 of [Nik79](@cite).

Note: for the implementation of the algorithm, we construct a lattice $T$ which
is unique in its genus, such that $D_T$ and $D_M(-1)$ are isometric and
$O(T)\to O(D_T)$ is surjective. We then classify all primitive
extensions of $M\oplus T$ modulo $O(D_T)$ (and modulo $O(M)$ for a
classification of primitive sublattices). To classify such primitive
extensions, we use Proposition 1.5.1 of [Nik79](@cite):

```@docs
primitive_extensions(::ZZLat, ::ZZLat)
```

We recall that a *primitive extension* of the orthogonal direct sum of two
integral integer lattices $M$ and $N$ is an overlattice $L$ of $M\oplus N$ such
that both $M$ and $N$ embed primitively in $L$ (via the natural embeddings
$M,N \to M\oplus N\subseteq L$). Such primitive extensions are obtained, and
classified, by classifying *gluings* between anti-isometric subgroups of the
respective discriminant groups of $M$ and $N$. The construction of an
overlattice is determined by the graph of such gluing.

## Equivariant primitive extensions

An *equivariant primitive extension* of a pair of integer lattices with
isometries $(M, f_M)$ and $(N, f_N)$ is a primitive extension of $M$ and $N$
obtained by gluing two subgroups which are respectively $D_{f_M}$ and
$D_{f_N}$ stable along a glue map which commutes with these two actions.
If such a gluing exists, then the overlattice $L$ of $M\oplus N$ is equipped with
an isometry $f_L$ which preserves both $M$ and $N$, and restricts to $f_M$ and
$f_N$ respectively.

```@docs
equivariant_primitive_extensions(::Union{ZZLatWithIsom, ZZLat}, ::Union{ZZLatWithIsom, ZZLat})
```

## Admissible equivariant primitive extensions

The following function is a major tool provided by [BH23](@cite). Given
a triple of even integer lattices with isometry $((A, a), (B, b), (C, c))$
and two prime numbers $p$ and $q$ (possibly equal), if $(A, B, C)$ is $p$-admissible,
this function returns representatives of isomorphism classes of equivariant primitive
extensions $(A, a)\oplus (B, b)\to (D, d)$ such that the type of $(D, d^q)$ is
equal to the type of $(C, c)$ (see [`type(::ZZLatWithIsom)`](@ref)).

```@docs
admissible_equivariant_primitive_extensions(::ZZLatWithIsom, ::ZZLatWithIsom, ::ZZLatWithIsom, ::IntegerUnion, ::IntegerUnion)
```
