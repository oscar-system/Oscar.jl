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

         $\mathcal R: K[V] \to K[V], f\to \frac{1}{|G|}\sum_{\pi\in G}(\pi f).$

!!! note
    We speak of *non-modular* invariant theory if $|G|$ is invertible in $K$, and of *modular* invariant theory otherwise.

!!! note
    In the non-modular case, using  Emmy Noether's result and the Reynolds operator, it is not too difficult to show that $K[V]^G$ is a free module over any of its graded Noether normalizations. That is, $K[V]^G$ is Cohen-Macaulay.

!!! note
    In the non-modular case, the Hilbert series of $K[V]^G$ can be precomputed via Molien's theorem. See [DK15](@cite) or [DJ98](@cite) for explicit formulas.

The algorithms for computing generators of invariant rings of finite groups proceed in two steps.
- First, compute a system of primary invariants.
- Then, compute a corresponding system of secondary invariants.

## Creating Invariant Rings

```@docs
invariant_ring(G::MatrixGroup)
```

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
