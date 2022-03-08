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

$(f \;\!   . \;\! \pi)  (x_1, \dots, x_n)  = f((x_1, \dots, x_n) \cdot \rho(\pi)) \text{ for all } \pi\in G.$


!!! note
    
    - The definition of linear reductivity guarantees the existence of a Reynolds operator $\mathcal R: K[V] \to K[V]$. 
    - By Hilbert's celebrated finiteness theorem, $K[V]^G$ is finitely generated as a $K$-algebra.
    - By a result of Hochster and Roberts, $K[V]^G$ is Cohen-Macaulay. 

In cases where the Reynold's operator can be explicitly handled, generators of invariant rings of linearly reductive groups can be found in two steps using Derksen's algorithm, see [Der99](@cite) :

- First, compute generators of Hilbert's null-cone ideal.
- Then, apply the Reynold's operator to these generators.

See also [DK15](@cite) and [DJ98](@cite).

## Creating Invariant Rings

There are no exact means to handle algebraically closed fields on the computer. For the computation of invariant rings in the above setting, on the other hand, there is no need to deal with explicit elements of ``G`` or with its group structure. The implementation of Derksen's algorithm in OSCAR can  handle situations where both $G$ and the representation $\rho$ are defined over an exact subfield $k$ of $K$ which is supported by OSCAR: 

- ``G`` is  specified as an affine algebraic variety by polynomials with coefficients in $k$;
- ``\rho: G \to \text{GL}(V) \cong \text{GL}_n(K)`` is specified by an $n\times n$ matrix whose entries are polynomials in the same variables as those specifying $G$, with coefficients in $k$.

All computations are then performed over $k$.

!!! warning
    OSCAR does neither check whether the affine variety defined by the given equations carries a group structure which makes it a linearly reductive group nor does it check whether the given $n\times n$ matrix really defines a representation.
    


## Basic Data Associated to Invariant Rings

## The Reynolds Operator

Omega-process

## Generators of Hilbert's Null-Cone Ideal

## Generators of the Invariant Ring

## Fundamental Systems of Invariants

## Invariant Rings as Affine Algebras






