```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["it_fg.md"]
```

# Invariants of Finite Groups

Using the notation from the introductory section, we now suppose that $\rho: G\to V$ is the representation
of a *finite* group $G$ on $V$. Recall that we identify $K[V]\cong K[x] $ via a fixed set of cooordinates $x_1, \dots, x_n\in V^*$.

!!! note
    - By Emmy Noether's finiteneness theorem, $K[V]^G$ is a finitely generated $K$-algebra of dimension $\dim K[V]^G = \dim K[V] = n$.
    - If the group order $|G|$ is invertible in $K$, then $K[V]^G$ is Cohen-Macaulay. In fact, in this case, $G$ is a linearly reductive group with explicitly given Reynolds operator

       $\mathcal R: K[V] \to K[V], f(x)\to \sum_{\pi\in G}(\pi f(x))$.

!!! note
    We speak of *non-modular* invariant theory if $|G|$ is invertible in $K$, and of *modular* invariant theory otherwise.

!!! note
    In the non-modular case, the Hilbert series of $K[V]^G$ is explicitly given as the Molien series of $G$, see [DK15](@cite) or [DJ98](@cite).

The algorithms for computing generators of invariant rings of finite groups proceed in two steps. First, compute a system of primary invariants. Then, compute a corresponding system of secondary invariants.
## Primary Invariants

## Secondary Invariants

## Invariant Rings and Fundamental Invariants
