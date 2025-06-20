# Schubert Calculus

To recall the definition of Schubert cycles on a Grassmannian $\mathrm{G}(k,n)$, we think of $\mathrm{G}(k, n)$ as the
Grassmannian $\mathrm{G}(k, W)$ of $k$-dimensional subspaces of an $n$-dimensional $K$-vector space $W$.

A *flag* in $W$ is a strictly increasing sequence of linear subspaces

$\{0\} \subset W_1 \subset \dots \subset W_{n-1} \subset W_n = W, \; \text{ with }\; \dim(W_i) = i.$

Let such a flag  $\mathcal{W}$ be given. For any sequence $a = (a_1, \ldots, a_k)$
of integers with $n-k \geq a_1 \geq \ldots \geq a_k \geq 0$, we define
the *Schubert cycle* $\Sigma_a(\mathcal{W})$ by setting

$\Sigma_a(\mathcal{W}):= \{ V \in \mathrm{G}(k, W)\mid \dim(W_{n-k+i-a_i} \cap V) \geq i, i = 1, \ldots, k \} \,.$

Then we have:
- The Schubert cycle $\Sigma_a(\mathcal{W})$ is an irreducible subvariety of $\mathrm{G}(k, W)$ of codimension $|a| = \sum a_i$.
- Its cycle class $[\Sigma_a(\mathcal{W})]$ does not depend on the choice of the flag.

We define the *Schubert class* of $a = (a_1, \ldots, a_k)$ to be the cycle class

$\sigma_a := [\Sigma_a(\mathcal{W})] \,.$

The number of Schubert classes on the Grassmannian $\mathrm{G}(k, n)$ is equal to $\binom{n}{k}$.

Instead of $\sigma_a$, we write $\sigma_{a_1,\ldots,a_s}$  whenever
$a = (a_1, \ldots, a_s, 0, \ldots, 0)$, and $\sigma_{p^i}$ whenever
$a = (p, \ldots, p, 0, \ldots, 0) \in \mathbb Z^i \times \{0\}^{k-i}$.

The classes  $\sigma_{1^i}$, $i = 1, \ldots, k$, and $\sigma_i$, $i = 1, \ldots, n-k$, are
called *special Schubert classes*. They are closely related to the Chern classes of the tautological
vector bundles on $\mathrm{G}(k, W)$. Recall:

- The *tautological subbundle* on $\mathrm{G}(k, W)$ is the vector bundle of rank $k$ whose fiber at $V \in \mathrm{G}(k, W)$ is the subspace $V \subset W$.
- The *tautological quotient bundle* on $\mathrm{G}(k, W)$ is the vector bundle of rank $(n-k)$ whose fiber at $V \in \mathrm{G}(k, W)$ is the quotient vector space $W/V$.

We denote these vector bundles by $S$ and $Q$, respectively.

The Chern classes of $S$ and $Q$ are

$c_i(S) = (-1)^i \sigma_{1^i} \; \text{ for } \;  i = 1, \ldots, k$

and

$c_i(Q) = \sigma_i \; \text{ for }\; i = 1, \ldots, n-k,$

respectively.

All Schubert classes form a $K$-vector space basis of the Chow ring of $\mathrm{G}(k,n)$.
The  Chern classes of $S$ (the special Schubert classes $\sigma_{1^i}$, $i=1, \ldots, k$)
form a minimal set of generators for the Chow ring of $\mathrm{G}(k,n)$ as a $K$-algebra.
See [EH16](@cite) for the relations on these generators.

```@docs
schubert_class(G::AbstractVariety, λ::Int...)
```

```@docs
schubert_classes(G::AbstractVariety)
```




