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

From a theoretical point of view, our basic setting in invariant theory consists of a group $G$ which acts on a vector space $V$ of finite
dimension $n$ over a field $K$. The action of $G$ on $V$ induces actions of $G$ on the dual vector space $V^\ast$,

$(\pi \cdot f)(v)=f(\pi^{-1}\cdot v),$

and, thus,  on the graded symmetric algebra

$K[V]:=S(V^*)=\bigoplus_{d\geq 0} S^d V^*.$

The fixed points of this action are the invariants of $G$, and the graded subalgebra

$K[V]^G=\{f\in K[V] \mid \pi\cdot  f=f {\text { for any }} \pi\in G\}\subset K[V]$

is its ring of invariants.

!!! note
    Except where mentioned otherwise, we will be in the favourable situation where $G$ is a linear reductive group which acts rationally on $V$. This has several important consequences:
    - There exists a Reynolds operator $\mathcal R: K[V] \rightarrow K[V]$. That is, $\mathcal R$ is a $K$-linear graded map which projects $K[V]$ onto $K[V]^G$, and which is a $K[V]^G$-module homomorphism;
    - by Hilbert's finiteness theorem, $K[V]^G$ is finitely generated as a $K$-algebra;
    - by a result of Hochster and Roberts, $K[V]^G$ is Cohen-Macaulay. Equivalently, $K[V]^G$ is a free module (of finite rank) over any of its Noether normalizations.
    If $k[V]^G$ is finitely generated as a $K$-algebra, we call any irredundant system of homogeneous generators a fundamental system of invariants of $k[V]^G$. By Nakayama's lemma, the number of elements in such a system is uniquely determined as the embedding dimension of $K[V]^G$. Similarly, the degrees of these elements are uniquely determined.


From a practical point of view, we will work with a fixed set of coordinates  $x_1, \dots, x_n\in V^*$, and $G$  will be a matrix group $G\subset \text{GL}_n(K) \cong \text{GL}_K(V^*) $, acting on $K[x_1, \dots, x_n]\cong K[V]$ by linear substitution:

$\pi\cdot f(x_1,\dots ,x_n) = f((x_1,\dots ,x_n)\cdot\pi).$

Accordingly, we will then write $K[x_1, \dots, x_n]^G\cong K[V]^G$.

The textbook

- [DK15](@cite)

and the survey article

- [DJ98](@cite)

provide details on theory and algorithms as well as references.
