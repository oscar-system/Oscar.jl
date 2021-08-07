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

Using the notation from the introductory section to this chapter, we now suppose that $G$ is a *finite* group.

!!! note
     - By a result of Emmy Noether, $K[V]$ is integral over $K[V]^G$. In particular,

          $\dim K[V]^G = \dim K[V] = n.$
         
         Moreover, $K[V]^G$ is finitely generated as a $K$-algebra.
		   
    - If the group order $|G|$ is invertible in $K$, then we have the explicit Reynolds operator

       $\mathcal R: K[V] \to K[V], f\to \frac{1}{|G|}\sum_{\pi\in G}(\pi \;\!  . f).$

!!! note
    We speak of *non-modular* invariant theory if $|G|$ is invertible in $K$, and of *modular* invariant theory otherwise.

!!! note
    In the non-modular case, using  Emmy Noether's result and the Reynolds operator, it is not too difficult to show that $K[V]^G$ is a free module over any of its graded Noether normalizations. That is, $K[V]^G$ is Cohen-Macaulay.

!!! note
    In the non-modular case, the Hilbert series of $K[V]^G$ can be precomputed via Molien's theorem. See [DK15](@cite) or [DJ98](@cite) for explicit formulas.

Having means to compute a $K$-basis for the invariants of each given degree, the algorithms for computing generators of invariant rings of finite groups proceed in two steps:

- First, compute a system of primary invariants $p_1,\dots, p_n$.
- Then, compute a system of secondary invariants with respect to $p_1,\dots, p_n$.

In the non-modular case, the Molien series allows one to precompute the number of $K$-linearly independent invariants for each given degree,

## Creating Invariant Rings

The invariant theory module of OSCAR  distinguishes two ways of how  finite groups and their actions on $K[x_1, \dots, x_n]\cong K[V]$ are given.

### Matrix Groups

Here, $G$ will be explicitly given as a matrix group $G\subset \text{GL}_n(K)\cong \text{GL}(V) $ by (finitely many) generating matrices, acting on $K[x_1, \dots, x_n]\cong K[V]$ by linear substitution:

$(\pi \;\!  .f) \;\! (x_1, \dots, x_n)  = f(\pi^{-1} \cdot (x_1, \dots, x_n)^T) \text{ for all } \pi\in G.$


```@docs
invariant_ring(G::MatrixGroup)
```

### Permutation Groups


## Basic Data Associated to Invariant Rings

## The Reynolds Operator

## Invariants of a Given Degree

## The Molien Series

## Primary Invariants

```@docs
primary_invariants(IR::InvRing)
```

## Secondary Invariants

```@docs
secondary_invariants(IR::InvRing)
```

## Fundamental Systems of Invariants

## Invariant Rings as Affine Algebras
