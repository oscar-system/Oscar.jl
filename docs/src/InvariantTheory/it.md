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

The basic setting in this chapter consists of a group $G$, a field $K$, a vector space
$V$ over $K$ of finite dimension $n,$  and a representation $\rho: G \to \text{GL}(V)$
of $G$ on $V$. Furthermore, $V^\ast$ denotes the dual vector space of $V$ and 
$x = \{x_1, \dots, x_n\}\subset V^*$ denotes a fixed set of cooordinates.

The action of $G$ on $V$ defines an action of  $G$ on $V^\ast$,

$(\pi f)(v)=f(\rho(\pi^{-1}) v).$

This extends naturally to an action of $G$ on the graded symmetric algebra

$K[x] = K[x_1, \dots, x_n] \cong K[V]:=S(V^*)=\bigoplus_{d\geq 0} S^d V^*.$

The *invariants* of $G$ are the fixed points of this action, its *ring of invariants* is the graded subalgebra

$K[x]^G = K[x_1, \dots, x_n]^G\cong  K[V]^G=\{f\in K[V] \mid \pi f=f {\text { for any }} \pi\in G\}\subset K[V].$

Clearly, $K[V]^G$ depends only on the image $\rho(G)\subset \text{GL}(V)$.

!!! note
    If $K[V]^G$ is finitely generated as a $K$-algebra, then any minimal system of homogeneous generators is called a *fundamental system of invariants* of $K[V]^G$. By Nakayama's lemma, the number of elements in such a system is uniquely determined as the embedding dimension of $K[V]^G$. Similarly, the degrees of these elements are uniquely determined.

!!! note
     If $K[V]^G$ is finitely generated as a $K$-algebra, then $K[V]^G$ admits a graded Noether normalization, that is, a Noether normalization $K[p_1, \dots, p_m] \subset K[V]^G$ with $p_1, \dots, p_m$ homogeneous. Given any such Noether normalization, $p_1, \dots, p_m$ is called a system of *primary invariants* of $K[V]^G$, and  any minimal system $s_0=1, s_1,\dots, s_l$ of homogeneous generators of $K[V]^G$ as a $K[p_1, \dots, p_m]$-module is called a system of *secondary invariants* of $K[V]^G$ with respect to $p_1, \dots, p_m$.

!!! note
    In all situations considered in this chapter, theoretical results will guarantee that $K[V]^G$ is finitely generated as a $K$-algebra. In addition, if not mentioned otherwise, the following will hold:
    - There exists a Reynolds operator $\mathcal R: K[V] \to K[V]$. That is, $\mathcal R$ is a $K$-linear graded map which projects $K[V]$ onto $K[V]^G$, and which is a $K[V]^G$-module homomorphism.
    - The ring $K[V]^G$ is Cohen-Macaulay. Equivalently, $K[V]^G$ is a free module (of finite rank) over any of its graded Noether normalizations.

The textbook

- [DK15](@cite)

and the survey article

- [DJ98](@cite)

provide details on theory and algorithms as well as references.
