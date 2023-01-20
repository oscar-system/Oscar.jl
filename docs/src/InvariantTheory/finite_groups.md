```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["finite_groups.md"]
```

# Invariants of Finite Groups

In this section, with notation as in the introduction to this chapter, $G$ will be a *finite* group.

!!! note
     The ssumption that $G$ is finite implies:
     - By a result of Emmy Noether, $K[V]$ is integral over $K[V]^G$. In particular,

          $\; \; \; \; \; \dim K[V]^G = \dim K[V] = n.$
         
         Moreover, $K[V]^G$ is finitely generated as a $K$-algebra.
		   
    - If the group order $|G|$ is invertible in $K$, then we have the explicit Reynolds operator

       $\; \; \; \; \; \mathcal R: K[V] \to K[V], f\mapsto \frac{1}{|G|}\sum_{\pi\in G}(f \;\!   . \;\! \pi).$

!!! note
    We speak of *non-modular* invariant theory if $|G|$ is invertible in $K$, and of *modular* invariant theory otherwise.

!!! note
    In the non-modular case, using  Emmy Noether's result and the Reynolds operator, it is not too difficult to show that $K[V]^G$ is a free module over any of its graded Noether normalizations. That is, $K[V]^G$ is Cohen-Macaulay. In the modular case, $K[V]^G$ may not be Cohen-Macaulay.

!!! note
    In the non-modular case, the Hilbert series of $K[V]^G$ can be precomputed as its Molien series. See [DK15](@cite) and [DJ98](@cite) for explicit formulas.

Knowing the Hilbert series means to know the dimension of each graded piece $K[V]^G_d$. This information can often be used to speed up algorithms for finding invariants.
The most basic task here is to compute the invariants of  some given degree $d$, that is, to find  an explicit $K$-basis of $K[V]^G_d$. There are two different approaches:

- The *Reynolds Operator Method*, available in  the non-modular case, applies the Reynolds operator to sufficiently many monomials in $K[x_1, \dots, x_n]_d\cong K[V]_d$,  and extracts a $K$-basis from the resulting generating set.
- The *Linear Algebra Method*, available in the non-modular and the modular case, finds the elements of a $K$-basis all at once by setting up and solving an appropriate $K$-linear system of equations.

These methods are, in particular, crucial to the computation of primary and secondary invariants. Primary invariants and irreducible secondary invariants together generate $K[V]^G$ as a $K$-algebra. Omitting redundant generators yields a system of fundamental invariants.
In the non-modular case, an alternative and typically more effective way to compute generators of $K[V]^G$ is King's algorithm which finds a system of fundamental invariants directly, without computing primary and secondary invariants. See [Kin13](@cite).

We discuss the relevant OSCAR functionality below.

## Creating Invariant Rings

### How Groups are Given

The invariant theory part of OSCAR  distinguishes two ways of how  finite groups and their actions on $K[x_1, \dots, x_n]\cong K[V]$ are specified:

#### Matrix Groups

Here, $G$ will be explicitly given as a matrix group $G\subset \text{GL}_n(K)\cong \text{GL}(V) $ by (finitely many) generating matrices, acting on $K[x_1, \dots, x_n]\cong K[V]$ by linear substitution:

$(f \;\!   . \;\! \pi)  (x_1, \dots, x_n)  = f((x_1, \dots, x_n) \cdot \rho(\pi)) \text{ for all } \pi\in G.$

#### Permutation Groups

Here, $G$ will be given as a permutation group, acting on $K[x_1, \dots, x_n]\cong K[V]$ by permuting the variables.

### Constructors for Invariant Rings

```@docs
invariant_ring(G::MatrixGroup)
```

## Basic Data Associated to Invariant Rings

If `IR` is the invariant ring $K[x_1,..., x_n]^G$ of a finite matrix group $G$, then

- `group(IR)` refers to $G$,
- `coefficient_ring(IR)` to $K$, and
- `polynomial_ring(IR)` to $K[x_1,..., x_n]$.

Moreover, `is_modular(IR)` returns `true` in the modular case, and
`false` otherwise.

###### Examples

```jldoctest
julia> K, a = CyclotomicField(3, "a")
(Cyclotomic field of order 3, a)

julia> M1 = matrix(K, [0 0 1; 1 0 0; 0 1 0])
[0   0   1]
[1   0   0]
[0   1   0]

julia> M2 = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1])
[1   0        0]
[0   a        0]
[0   0   -a - 1]

julia> G = MatrixGroup(3, K, [ M1, M2 ])
Matrix group of degree 3 over Cyclotomic field of order 3

julia> IR = invariant_ring(G)
Invariant ring of
Matrix group of degree 3 over Cyclotomic field of order 3
with generators
AbstractAlgebra.Generic.MatSpaceElem{nf_elem}[[0 0 1; 1 0 0; 0 1 0], [1 0 0; 0 a 0; 0 0 -a-1]]

julia> group(IR)
Matrix group of degree 3 over Cyclotomic field of order 3

julia> coefficient_ring(IR)
Cyclotomic field of order 3

julia> R = polynomial_ring(IR)
Multivariate Polynomial Ring in x[1], x[2], x[3] over Cyclotomic field of order 3 graded by
  x[1] -> [1]
  x[2] -> [1]
  x[3] -> [1]

julia> x = gens(R)
3-element Vector{MPolyElem_dec{nf_elem, AbstractAlgebra.Generic.MPoly{nf_elem}}}:
 x[1]
 x[2]
 x[3]

julia> is_modular(IR)
false

```

## The Reynolds Operator

```@docs
reynolds_operator(IR::InvRing{FldT, GrpT, T}, f::T) where {FldT, GrpT, T <: MPolyElem}

reynolds_operator(IR::InvRing{FldT, GrpT, T}, f::T, chi::GAPGroupClassFunction) where {FldT, GrpT, T <: MPolyElem}
```

## Invariants of a Given Degree

```@docs
basis(IR::InvRing, d::Int, algo::Symbol = :default)

basis(IR::InvRing, d::Int, chi::GAPGroupClassFunction)
```

```@docs
iterate_basis(IR::InvRing, d::Int, algo::Symbol = :default)

iterate_basis(IR::InvRing, d::Int, chi::GAPGroupClassFunction)
```

## The Molien Series

```@docs
 molien_series([S::PolyRing], I::InvRing, [chi::GAPGroupClassFunction])
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
fundamental_invariants(IR::InvRing, algo::Symbol = :default; beta::Int = 0)
```

## Invariant Rings as Affine Algebras

```@docs
affine_algebra(IR::InvRing)
```

## Semi-invariants / relative invariants
```@docs
semi_invariants(IR::InvRing, chi::GAPGroupClassFunction)
```
