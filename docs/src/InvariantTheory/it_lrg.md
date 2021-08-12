```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["it_lrg.md"]
```

# Invariants of Linearly Reductive Groups

In this section, with notation as in the introduction to this chapter, $G$ will be a *linearly algebraic group* over an algebraically closed field $K$, and ``\rho: G \to \text{GL}(V)\cong \text{GL}_n(K)`` will be a *rational* representation of $G$. As in the previous sections, ``G`` will act on $K[V]\cong K[x_1, \dots, x_n]$ by linear substitution:

$(\pi \;\!  . \;\! f) \;\! (x_1, \dots, x_n)  = f(\rho(\pi^{-1}) \cdot (x_1, \dots, x_n)^T) \text{ for all } \pi\in G.$

!!! note
    
    - By the very definition of linear reductivity, there is a Reynolds operator $\mathcal R: K[V] \to K[V]$. 
    - By Hilbert's celebrated finiteness theorem, $K[V]^G$ is finitely generated as a $K$-algebra.
    - By a result of Hochster and Roberts, $K[V]^G$ is Cohen-Macaulay. 

Provided that a Reynold's operator is explicitly given, generators of invariant rings of linearly reductive groups can be found in two steps using Derksen's algorithm:

- First, compute generators of Hilbert's null-cone ideal.
- Then, apply the Reynold's operator to these generators.

See [DK15](@cite) and [DJ98](@cite).

## Creating Invariant Rings

The linearly algebraic groups dealt with by the invariant theory part of OSCAR  are defined over an exact subfield $k$ of $K$ which is supported by OSCAR. That is:

- ``G`` is  specified as an algebraic subgroup of $\text{GL}_t(K)$ by polynomial equations over $k$,  for some $t$.
- ``\rho: G \to \text{GL}(V) \cong \text{GL}_n(K)`` is a *rational* representation of $G$ given by polynomials over $k$.

!!! note
    There are no exact means to handle algebraically closed fields on the computer. For the computation of invariant rings in the above setting, however, there is no need to deal with explicit elements of ``G``. All computations are performed over $k$.

## Basic Data Associated to Invariant Rings

## The Reynolds Operator

Omega-process

## Generators of Hilbert's Null-Cone Ideal

## Generators of the Invariant Ring

## Fundamental Systems of Invariants

## Invariant Rings as Affine Algebras






