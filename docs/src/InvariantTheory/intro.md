```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["intro.md"]
```

# Introduction

The invariant theory part of OSCAR provides functionality for computing polynomial invariants
of group actions, focusing on finite and linearly reductive groups, respectively.

The basic setting in this context consists of a group $G$, a field $K$, a vector space
$V$ over $K$ of finite dimension $n,$ and  a representation $\rho: G \to \text{GL}(V)$ of $G$ on $V$.
The induced action on the dual vector space $V^\ast$,

$V^\ast  \times G \rightarrow V^\ast, (f, \pi)\mapsto f \;\!   . \;\! \pi  := f\circ \rho(\pi),$
 
extends to an action of $G$ on the graded symmetric algebra

$K[V]:=S(V^*)=\bigoplus_{d\geq 0} S^d V^*$

which preserves the grading.

The *invariants* of $G$ are the fixed points of this action, its *invariant ring* is the graded subalgebra

$K[V]^G:=\{f\in K[V] \mid f \;\!   . \;\! \pi =f {\text { for any }} \pi\in G\} \subset K[V].$
       	   
Explicitly, the choice of a basis of $V$ and its dual basis, say, $\{x_1, \dots, x_n\}$ of $V^*$
gives rise to isomorphisms $\text{GL}(V) \cong \text{GL}_n(K)$ and $K[V]\cong  K[x_1, \dots, x_n]$.
After identifying $\text{GL}(V)$ with $\text{GL}_n(K)$ and $K[V]$ with $K[x_1, \dots, x_n]$ by means of
these isomorphisms, the action of $G$ on $K[V]$ is given as follows:

$(f \;\!   . \;\! \pi)  (x_1, \dots, x_n)  = f((x_1, \dots, x_n) \cdot \rho(\pi)).$

Accordingly, $K[V]^G$ may be regarded as a graded subalgebra of $K[x_1, \dots, x_n]$:

$K[V]^G \cong K[x_1, \dots, x_n]^G :=\{f\in K[x_1, \dots, x_n] \mid f \;\!   . \;\! \pi =f {\text { for any }} \pi\in G\}.$

The main objective of invariant theory in OSCAR is the computation of $K$-algebra generators for invariant rings.

!!! note
    If $K[V]^G$ is finitely generated as a $K$-algebra, then any minimal system of homogeneous generators is called a *fundamental system of invariants* for $K[V]^G$. By Nakayama's lemma, the number of elements in such a system is uniquely determined as the embedding dimension of $K[V]^G$. Similarly, the degrees of these elements are uniquely determined.

!!! note
    If $K[V]^G$ is finitely generated as a $K$-algebra, then $K[V]^G$ admits a graded Noether normalization, that is, a Noether normalization $K[p_1, \dots, p_m] \subset K[V]^G$ with $p_1, \dots, p_m$ homogeneous. Given any such Noether normalization, $p_1, \dots, p_m$ is called a *homogeneous system of parameters* or a *system of primary invariants* for $K[V]^G$, and  any minimal system $s_0=1, s_1,\dots, s_l$ of homogeneous generators of $K[V]^G$ as a $K[p_1, \dots, p_m]$-module is called a *system of secondary invariants* for $K[V]^G$ with respect to $p_1, \dots, p_m$. A secondary invariant $s_i\neq 1$ is called *irreducible* if it cannot be written as a polynomial expression in the primary invariants and the other secondary invariants. The  irreducible secondary invariants form a minimal system of homogeneous generators for $K[V]^G$ as a $K[p_1, \dots, p_m]$-algebra. Somewhat abusing notation, we call every minimal system of homogeneous generators for $K[V]^G$ as a $K[p_1, \dots, p_m]$-algebra a *system of irreducible secondary invariants*.


!!! note
    For the invariant rings handled by OSCAR, the assumption that $K[V]^G$ is finitely generated as a $K$-algebra will be guaranteed by theoretical results. In addition, where not mentioned otherwise, the following will hold:
    - There exists a Reynolds operator $\mathcal R: K[V] \to K[V]$. That is, $\mathcal R$ is a $K$-linear graded map which projects $K[V]$ onto $K[V]^G$, and which is a $K[V]^G$-module homomorphism.
    - The ring $K[V]^G$ is Cohen-Macaulay. Equivalently, $K[V]^G$ is a free module (of finite rank) over any of its graded Noether normalizations.

The textbook

- [DK15](@cite)

and the survey article

- [DJ98](@cite)

provide details on theory and algorithms as well as references.
