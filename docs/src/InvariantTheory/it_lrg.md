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

In this section, with notation as in the introduction to this chapter, $G$ will always be a *linearly algebraic group* over an algebraically closed field $K$, and ``\rho: G \to \text{GL}(V)`` will be a *rational* representation of $G$. As in the previous sections, ``G`` will act on $K[x_1, \dots, x_n]\cong K[V]$ by linear substitution:

$(\pi \;\!  .f) \;\! (x_1, \dots, x_n)  = f(\rho(\pi^{-1}) \cdot (x_1, \dots, x_n)^T) \text{ for all } \pi\in G.$

!!! note
    
    - By the very definition of linear reductivity, there is a Reynolds operator $\mathcal R: K[V] \to K[V]$. 
    - By Hilbert's celebrated finiteness theorem, $K[V]^G$ is finitely generated as a $K$-algebra.
    - By a result of Hochster and Roberts, $K[V]^G$ is Cohen-Macaulay. 

## Creating Invariant Rings

The invariant theory module of OSCAR handles linearly algebraic groups which are defined over an exact subfield $k$ of $K$ which is supported by OSCAR. That is:

- ``G`` is  specified as an algebraic subgroup of $\text{GL}_t(K)$ by polynomial equations over $k$,  for some $t$.
- ``\rho: G \to \text{GL}_n(K)\cong \text{GL}(V)`` is a *rational* representation of $G$ given by polynomials over $k$.

!!! note
    There are no exact means to handle algebraically closed fields on the computer. In the above setting, we do not deal with explicit elements of ``G``.


## Basic Data Associated to Invariant Rings

## The Reynolds Operator

Omega-process

## Generators of the Null-Cone

## Generators of the Invariant Ring

## Fundamental Systems of Invariants

## Invariant Rings as Affine Algebras






