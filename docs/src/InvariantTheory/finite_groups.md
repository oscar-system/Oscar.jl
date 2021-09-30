```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["finite_groups.md"]
```

# Invariants of Finite Groups

In this section, with notation as in the introduction to this chapter, $G$ will always be a *finite* group.

!!! note
     - By a result of Emmy Noether, $K[V]$ is integral over $K[V]^G$. In particular,

          $\dim K[V]^G = \dim K[V] = n.$
         
         Moreover, $K[V]^G$ is finitely generated as a $K$-algebra.
		   
    - If the group order $|G|$ is invertible in $K$, then we have the explicit Reynolds operator

       $\mathcal R: K[V] \to K[V], f\mapsto \frac{1}{|G|}\sum_{\pi\in G}(\pi \;\!  . \;\! f).$

!!! note
    We speak of *non-modular* invariant theory if $|G|$ is invertible in $K$, and of *modular* invariant theory otherwise.

!!! note
    In the non-modular case, using  Emmy Noether's result and the Reynolds operator, it is not too difficult to show that $K[V]^G$ is a free module over any of its graded Noether normalizations. That is, $K[V]^G$ is Cohen-Macaulay.

!!! note
    In the non-modular case, the Hilbert series of $K[V]^G$ can be precomputed as its Molien series. See [DK15](@cite) and [DJ98](@cite) for explicit formulas.

Having means to compute a $K$-basis for the invariants of each given degree, the algorithms for computing generators of invariant rings of finite groups proceed in two steps:

- First, compute a system of primary invariants $p_1,\dots, p_n$.
- Then, compute a system of secondary invariants with respect to $p_1,\dots, p_n$.

In the non-modular case, the Molien series allows one to precompute the number of $K$-linearly independent invariants for each given degree,

## Creating Invariant Rings

The invariant theory part of OSCAR  distinguishes two ways of how  finite groups and their actions on $K[x_1, \dots, x_n]\cong K[V]$ are specified.

### Matrix Groups

Here, $G$ will be explicitly given as a matrix group $G\subset \text{GL}_n(K)\cong \text{GL}(V) $ by (finitely many) generating matrices, acting on $K[x_1, \dots, x_n]\cong K[V]$ by linear substitution:

$(\pi \;\!  . \;\! f) \;\! (x_1, \dots, x_n)  = f(\pi^{-1} \cdot (x_1, \dots, x_n)^T) \text{ for all } \pi\in G.$


```@docs
invariant_ring(G::MatrixGroup)
```

### Permutation Groups


## Basic Data Associated to Invariant Rings

If `IR` is the invariant ring $K[x_1,..., x_n]^G$ of a finite matrix group $G$, then

- `group(IR)` refers to $G$,
- `coefficient_ring(IR)` to $K$, and
- `polynomial_ring(IR)` to $K[x_1,..., x_n]$.

Moreover, `ismodular(IR)` returns `true` in the modular case, and
`false` otherwise.

###### Examples

```@repl oscar
K, a = CyclotomicField(3, "a")
M1 = matrix(K, [0 0 1; 1 0 0; 0 1 0])
M2 = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1])
G = MatrixGroup(3, K, [ M1, M2 ])
IR = invariant_ring(G)
group(IR)
coefficient_ring(IR)
R = polynomial_ring(IR)
x=gens(R)
ismodular(IR)
```

## The Reynolds Operator

```@docs
reynolds_operator(IR::InvRing{FldT, GrpT, T}, f::T) where {FldT, GrpT, T <: MPolyElem}
```

## Invariants of a Given Degree

```@docs
basis(IR::InvRing, d::Int)
```

## The Molien Series

```@docs
molien_series(IR::InvRing)
```

## Primary Invariants

```@docs
primary_invariants(IR::InvRing)
```

## Secondary Invariants

```@docs
secondary_invariants(IR::InvRing)
```

```@docs
irreducible_secondary_invariants(IR::InvRing)
```

## Fundamental Systems of Invariants

```@docs
fundamental_invariants(IR::InvRing)
```

## Invariant Rings as Affine Algebras

```@docs
affine_algebra(IR::InvRing)
```
