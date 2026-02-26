```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Schubert calculus

The term *Schubert calculus* refers to the intersection theory of Grassmannians. This is to honor Hermann Schubert,
a German 18th century mathematician who solved a variety of problems in enumerative geometry by reducing them
to the combinatorics and intersection theory of certain cycles on Grassmannians. These cycles are nowadays called
*Schubert cycles*, their cycle classes are called *Schubert classes*.

We show how to define such classes in OSCAR.

To recall the definition of Schubert cycles on a Grassmannian $\mathrm{Gr}(k,n)$, we think of $\mathrm{Gr}(k, n)$ as the
Grassmannian $\mathrm{Gr}(k, W)$ of $k$-dimensional subspaces of an $n$-dimensional $K$-vector space $W$.

A *flag* in $W$ is a strictly increasing sequence of linear subspaces

$\{0\} \subset W_1 \subset \dots \subset W_{n-1} \subset W_n = W, \; \text{ with }\; \dim(W_i) = i.$

Let such a flag $\mathcal{W}$ be given. For any sequence $a = (a_1, \ldots, a_k)$
of integers with $n-k \geq a_1 \geq \ldots \geq a_k \geq 0$, we define
the *Schubert cycle* $\Sigma_a(\mathcal{W})$ by setting

$\Sigma_a(\mathcal{W}):= \{ V \in \mathrm{Gr}(k, W)\mid \dim(W_{n-k+i-a_i} \cap V) \geq i, i = 1, \ldots, k \} \,.$

Then we have:
- The Schubert cycle $\Sigma_a(\mathcal{W})$ is an irreducible subvariety of $\mathrm{Gr}(k, W)$ of codimension $|a| = \sum a_i$.
- Its cycle class $[\Sigma_a(\mathcal{W})]$ does not depend on the choice of the flag.

We define the *Schubert class* of $a = (a_1, \ldots, a_k)$ to be the cycle class

$\sigma_a := [\Sigma_a(\mathcal{W})] \,.$

The number of Schubert classes on the Grassmannian $\mathrm{Gr}(k, n)$ is equal to $\binom{n}{k}$.

Instead of $\sigma_a$, we write $\sigma_{a_1,\ldots,a_s}$ whenever
$a = (a_1, \ldots, a_s, 0, \ldots, 0)$, and $\sigma_{p^i}$ whenever
$a = (p, \ldots, p, 0, \ldots, 0) \in \mathbb Z^i \times \{0\}^{k-i}$.

The classes $\sigma_{1^i}$, $i = 1, \ldots, k$, and $\sigma_i$, $i = 1, \ldots, n-k$, are
called *special Schubert classes*. They are closely related to the Chern classes of the tautological
vector bundles on $\mathrm{Gr}(k, W)$. Recall:

- The *tautological subbundle* on $\mathrm{Gr}(k, W)$ is the vector bundle of rank $k$ whose fiber at $V \in \mathrm{Gr}(k, W)$ is the subspace $V \subset W$.
- The *tautological quotient bundle* on $\mathrm{Gr}(k, W)$ is the vector bundle of rank $(n-k)$ whose fiber at $V \in \mathrm{Gr}(k, W)$ is the quotient vector space $W/V$.

We denote these vector bundles by $S$ and $Q$, respectively.

The Chern classes of $S$ and $Q$ are

$\operatorname{c}_i(S) = (-1)^i \sigma_{1^i} \; \text{ for } \; i = 1, \ldots, k$

and

$\operatorname{c}_i(Q) = \sigma_i \; \text{ for }\; i = 1, \ldots, n-k,$

respectively.

All Schubert classes form a $\mathbb{Q}$-vector space basis of the Chow ring of $\mathrm{Gr}(k,n)$.
The Chern classes of $S$ (the special Schubert classes $\sigma_{1^i}$, $i=1, \ldots, k$)
form a minimal set of generators for the Chow ring of $\mathrm{Gr}(k,n)$ as a $K$-algebra.
See [EH16](@cite) for the relations on these generators.

## Some classic computations

We illustrate Schubert calculus with a few classic enumerative problems.

**How many lines in $\mathbb P^3$ meet four general lines?**

The Grassmannian of lines in $\mathbb P^3$ is $\mathrm{Gr}(2,4)$.
The Schubert class $\sigma_1$ represents lines meeting a given line.
The answer is obtained by integrating $\sigma_1^4$:

```jldoctest
julia> G = abstract_grassmannian(2, 4)
AbstractVariety of dim 4

julia> s1 = schubert_class(G, 1)
-c[1]

julia> integral(s1^4)
2

```

**How many lines lie on a cubic surface in $\mathbb P^3$?**

To count lines on a cubic surface in $\mathbb P^3$,
we consider the zero locus of a section of $\operatorname{Sym}^3 S^*$ on $\mathrm{Gr}(2,4)$.
Here $S$ is the tautological subbundle, and a line lies on the surface
if and only if the section vanishes on it.
The answer is the top Chern class of $\operatorname{Sym}^3 S^*$:

```jldoctest
julia> S, Q = tautological_bundles(G);

julia> integral(top_chern_class(symmetric_power(dual(S), 3)))
27

```

**How many lines lie on a quintic threefold in $\mathbb P^4$?**

The Grassmannian of lines in $\mathbb{P}^4$ is $\mathrm{Gr}(2,5)$.
To count lines on the general quintic threefold,
we consider the zero locus of a section of $\operatorname{Sym}^5 S^*$:

```jldoctest
julia> G = abstract_grassmannian(2, 5)
AbstractVariety of dim 6

julia> S = tautological_bundles(G)[1];

julia> integral(top_chern_class(symmetric_power(dual(S), 5)))
2875

```

## Functions

```@docs
schubert_class(G::AbstractVariety, λ::Int...)
```

```@docs
schubert_classes(G::AbstractVariety)
```
