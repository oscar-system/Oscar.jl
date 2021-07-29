```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["it.md"]
```

# Introduction

Our basic setting in invariant theory consists of a group $G$, a field $K$, 
a vector space $V$ over $K$ of finite dimension $n,$  and a representation $\rho: G \to \text{GL}(V)$
of $G$ on $V$. We write $V^\ast$ for the dual vector space of $V$ and suppose that a fixed set of
cooordinates $x = \{x_1, \dots, x_n\}\in V^*$ is chosen.


The action of $G$ on $V$ induces an action of  $G$ on $V^\ast$: Set

$(\pi f)(v)=f(\rho(\pi^{-1}) v).$

This, in turn, defines an action of $G$ on the graded symmetric algebra

$K[x] = K[x_1, \dots, x_n] \cong  K[V]=S(V^*)=\bigoplus_{d\geq 0} S^d V^*.$

The *invariants* of $G$ are the fixed points of this action, its *ring of invariants* is the graded subalgebra

$K[x]^G\cong K[V]^G=\{f\in K[V] \mid \pi f=f {\text { for any }} \pi\in G\}\subset K[V].$

Clearly, this ring depends only on the image $\rho(G)\subset \text{GL}(V)$.

!!! note
    If not mentioned otherwise, we will be in the favourable situation where $G$ is a linear reductive group which acts rationally on $V$. This has several important consequences:
    - There exists a Reynolds operator $\mathcal R: K[V] \to K[V]$. That is, $\mathcal R$ is a $K$-linear graded map which projects $K[V]$ onto $K[V]^G$, and which is a $K[V]^G$-module homomorphism.
    - By Hilbert's finiteness theorem, $K[V]^G$ is finitely generated as a $K$-algebra.
    - By a result of Hochster and Roberts, $K[V]^G$ is Cohen-Macaulay. Equivalently, $K[V]^G$ is a free module (of finite rank) over any of its Noether normalizations.

!!! note
    If $k[V]^G$ is finitely generated as a $K$-algebra, then any minimal system of homogeneous generators is called a *fundamental system of invariants* of $k[V]^G$. By Nakayama's lemma, the number of elements in such a system is uniquely determined as the embedding dimension of $K[V]^G$. Similarly, the degrees of these elements are uniquely determined.

!!! note
    If $K[V]^G$ is finitely generated as a $K$-algebra, and $K[p_1, \dots, p_m] \subset K[V]^G$ is any Noether normalization, then $p_1, \dots, p_m$ is called a system of *primary invariants*. Given such a system $p_1, \dots, p_m$, we call any minimal system $s_0=1, s_1,\dots, s_l$ of homogeneous generators of $K[V]^G$ as a $K[p_1, \dots, p_m]$-module a system of *secondary invariants*.


The textbook

- [DK15](@cite)

and the survey article

- [DJ98](@cite)

provide details on theory and algorithms as well as references.
