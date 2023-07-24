```@meta
CurrentModule = Oscar
```

We introduce here the necessary definitions and results which lie behind the
methods presented. Most of the content is taken from [Nik79](@cite).

# Primitive embeddings in even lattices

## Nikulin's theory

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
into an even lattice in a given genus (see Proposition 1.15.1 in [Nik79](@cite)).
In the latter proposition, it is explained how to classify such embeddings as
isomorphic embeddings or as isomorphic sublattices.

Such a method can be algorithmically implemented, however it tends to be slow
and inefficient in general for large rank or determinant. But, in the case
where the discriminant groups are (elementary) $p$-groups, the method can be
more efficient.

The ultimate goal of the project is to make all kind of computations of
primitive embeddings available. For now, we only cover the case where one of
the two lattices involved is *$p$-primary*, i.e. its discriminant is an abelian
$p$-group. Note that this covers the case of unimodular lattices, of course.
We provide 4 kinds of output:
* A boolean, which only returns whether there exists a primitive embedding;
* A single primitive embedding as soon as the algorithm computes one;
* A list of representatives of isomorphism classes of primitive embeddings;
* A list of representatives of isomorphism classes of primitive sublattices.

```@docs
primitive_embeddings_in_primary_lattice(::ZZLat, ::ZZLat)
primitive_embeddings_of_primary_lattice(::ZZLat, ::ZZLat)
```

Note that the previous two functions require the first lattice of the input to be
unique in its genus. Otherwise, one can specify a genus, or its invariants, as a
first input:

```@docs
primitive_embeddings_in_primary_lattice(::ZZGenus, ::ZZLat)
primitive_embeddings_in_primary_lattice(::TorQuadModule, ::Tuple{Int, Int},
::ZZLat)
primitive_embeddings_of_primary_lattice(::ZZGenus, ::ZZLat)
primitive_embeddings_of_primary_lattice(::TorQuadModule, ::Tuple{Int, Int},
::ZZLat)
```

In order to compute such primitive embeddings of a lattice `M` into a lattice
`L`, one first computes the possible genera for the orthogonal of `M` in `L`
(after embedding), and for each lattice `N` in such a genus, one computes
isomorphism classes of *primitive extensions* of $M \perp N$ modulo $\bar{O}(N)$
(and $\bar{O}(M)$ in the case of classification of primitive sublattices of `L`
isometric to `M`).

We recall that a *primitive extension* of the orthogonal direct sum of two
integral integer lattices `M` and `N` is an overlattice `L` of $M\perp N$ such
that both `M` and `N` embed primitively in `L` (via the natural embeddings
$M,N \to M\perp N\subseteq L$). Such primitive extensions are obtained, and
classified, by looking for *gluings* between anti-isometric subgroups of the
respective discriminant groups of `M` and `N`. The construction of an
overlattice is determined by the graph of such glue map. 

## Admissible equivariant primitive extensions

The following function is an interesting tool provided by [BH23](@cite). Given
a triple of integer lattices with isometry `((A, a), (B, b), (C, c))` and two prime
numbers `p` and `q` (possibly equal), if `(A, B, C)` is `p`-admissible, this
function returns representatives of isomorphism classes of equivariant primitive
extensions $(A, a)\perp (B, b)\to (D, d)$ such that the type of $(D, d^p)$ is
equal to the type of $(C, c)$ (see [`type(::ZZLatWithIsom)`](@ref)).

```@docs
admissible_equivariant_primitive_extensions(::ZZLatWithIsom, ::ZZLatWithIsom, ::ZZLatWithIsom, ::Integer, ::Integer)
```

An *equivariant primitive extension* of a pair of integer lattices with
isometries $(M, f_M)$ and $(N, f_N)$ is a primitive extension of `M` and `N`
obtained by gluing two subgroups which are respectively $\bar{f_M}$ and
$\bar{f_N}$ stable along a glue map which commutes with these two actions. If
such a gluing exists, then the overlattice `L` of $M\perp N$ is equipped with
an isometry $f_L$ which preserves both `M` and `N`, and restricts to $f_M$ and
$f_N$ respectively.
