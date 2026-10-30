# [Introduction](@id invariant_theory)

The invariant theory part of OSCAR provides functionality for computing polynomial invariants
of group actions, focusing on finite groups, linearly reductive groups, and tori, respectively.

The basic setting in this context consists of a group $G$, a field $K$, a vector space
$V$ over $K$ of finite dimension $n,$ and  a representation $\rho: G \to \text{GL}(V)$ of $G$ on $V$.
The induced right action on the dual vector space $V^\ast$,

$V^\ast  \times G \to V^\ast, (f, \pi)\mapsto f \;\!   . \;\! \pi  := f\circ \rho(\pi),$

extends to a right action of $G$ on the graded symmetric algebra

$K[V]:=S(V^*)=\bigoplus_{d\geq 0} S^d V^*$

which preserves the grading.

!!! note
    In OSCAR, group actions are by convention assumed to be right actions and we follow this convention with our definition above.
    Note, however, that the left action given by $\pi \;\!   . \;\! f :=  f \circ \rho(\pi^{-1})$ is quite common in the literature.
	
The *invariants* of $G$ are the fixed points of the action defined above, its *invariant ring* is the graded subalgebra

$K[V]^G:=\{f\in K[V] \mid f \;\!   . \;\! \pi =f {\text { for any }} \pi\in G\} \subset K[V].$

Explicitly, fixing a basis of $V$ and its dual basis, say, $\{x_1, \dots, x_n\}$ of $V^*$,
we may identify $\operatorname{GL}(V) \cong \operatorname{GL}_n(K)$ and $K[V]\cong K[x_1, \dots, x_n]$.
Then the action of an element $\pi \in G$ with $\rho(\pi) = (a_{i, j})$ on a polynomial $f\in K[x_1,\dots, x_n]$ is given as follows:

$(f \;\!   . \;\! \pi)  (x_1, \dots, x_n) = f\bigl(\sum_j a_{1, j}x_j, \dots, \sum_j a_{n, j}x_j\bigr).$

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


## Comparison with other sources

The left action $\pi \;\!   . \;\! f :=  f \circ \rho(\pi^{-1})$ has the same invariants as the right action defined above.
The invariants do change if the transpose of $\rho(\pi)$ is substituted instead, that is, $x_j \mapsto \sum_i a_{i, j}x_i$: this yields the invariants of the dual representation.

- [Stu93](@cite), Section 1.3, maps $f$ to $f \circ \pi$, as OSCAR does.
- Magma (`InvariantRing`, `f^a`) and Singular (`invariant_ring` from `finvar.lib`) substitute $\rho(\pi)$ as OSCAR does.
- [DK15](@cite), Section 2.1, defines $(\pi \cdot f)(v) := f(\pi^{-1} \cdot v)$, which agrees with OSCAR if matrices act on column vectors. With explicit matrices, however, it substitutes the transpose: the columns of a matrix hold the images of the $x_i$ (Example 3.9.6), and the graph ideal of Algorithm 4.1.9 uses $y = A^t x$ (Section 4.1.3).

For example, OSCAR finds the invariant $x_2$ for the group generated by `[1 1; 0 1]` over $\mathbb{F}_3$, whereas [DK15](@cite), Example 3.9.13 with $q = 3$, gives $x_1$, which OSCAR finds for the transposed matrix:

```jldoctest
julia> M = matrix(GF(3), [1 1; 0 1]);

julia> fundamental_invariants(invariant_ring(matrix_group(M)))
2-element Vector{MPolyDecRingElem{FqFieldElem, FqMPolyRingElem}}:
 x[2]
 x[1]^3 + 2*x[1]*x[2]^2

julia> fundamental_invariants(invariant_ring(matrix_group(transpose(M))))
2-element Vector{MPolyDecRingElem{FqFieldElem, FqMPolyRingElem}}:
 x[1]
 x[1]^2*x[2] + 2*x[2]^3
```


## Contact

Please direct questions about this part of OSCAR to the following people:
* [Wolfram Decker](https://math.rptu.de/en/wgs/agag/people/head/decker),
* [Max Horn](https://math.rptu.de/en/wgs/agag/people/head/prof-dr-max-horn),
* [Johannes Schmitt](https://joschmitt.eu/).

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).
