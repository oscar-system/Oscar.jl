# Invariants of Linearly Reductive Groups

In this section, with notation as in the [introduction](@ref invariant_theory) to this chapter,
$G$ will be a *linear algebraic group* over an algebraically closed
field $K$, $\rho: G \to \text{GL}(V)\cong \text{GL}_n(K)$ will
be a *rational* representation of $G$, and
$G$ will act on $K[V]\cong K[x_1, \dots, x_n]$ by linear
substitution: If $\rho(\pi) = (a_{i, j})$, then 

$(f \;\!   . \;\! \pi)  (x_1, \dots, x_n) = f\bigl(\sum_j a_{1, j}x_j, \dots, \sum_j a_{n, j}x_j\bigr).$

!!! note
    
    - The definition of linear reductivity guarantees the existence of a Reynolds operator $\mathcal R: K[V] \to K[V]$. 
    - By Hilbert's celebrated finiteness theorem, $K[V]^G$ is finitely generated as a $K$-algebra.
    - By a result of Hochster and Roberts, $K[V]^G$ is Cohen-Macaulay.

In cases where the Reynolds operator is explicitly known, generators of invariant rings of linearly reductive groups
can be found in two steps using Derksen's algorithm, see [Der99](@cite) :

- First, compute generators of Hilbert's null-cone ideal.
- Then, apply the Reynolds operator to these generators.

See also [DK15](@cite) and [DJ98](@cite).

## Creating Invariant Rings

### How Linearly Reductive Groups and Their Representations are Given

For the computation of invariant rings in the above setting, there is no need to deal with explicit elements of ``G`` or with its group structure.
The implementation of Derksen's algorithm in OSCAR can  handle situations where both $G$ and the representation $\rho$ are defined over an exact
subfield $k$ of $K$ which is supported by OSCAR: 

- ``G`` is  specified as an affine algebraic variety by polynomials with coefficients in $k$;
- ``\rho: G \to \text{GL}(V) \cong \text{GL}_n(K)`` is specified by an $n\times n$ matrix whose entries are polynomials in the same variables as those specifying $G$, with coefficients in $k$.

!!! note
    Proceeding as above is not a problem: Derksen's algorithm relies on Gröbner bases techniques and means to compute
    Reynolds operators. It does, thus, not change the initial ground field $k$. That is, all computations are performed over $k$
	and computations over any extension field of $k$ would lead to the same results.

In OSCAR, the basic set-up for a linearly reductive group in the context of Derksen's algorithm is provided by the 
function `linearly_reductive_group`. At current state, this only supports rational actions of  the special linear group
$\text{SL}_m$ in characteristic zero (for $m \geq 2$, $\text{SL}_m$ is not linearly reductive in positive
characteristic). For the action of this group by linear
substitution on $m$-ary forms of degree $d$, an explicit Reynolds operator is
given by Cayley's Omega-process. We show this at work later in this section.


```@docs
linearly_reductive_group(sym::Symbol, m::Int, k::Field)
```

```@docs
linearly_reductive_group(sym::Symbol, m::Int, R::MPolyRing)
```

```@docs
linearly_reductive_group(I::MPolyIdeal)
```

```@docs
representation_reductive_group(G::LinearlyReductiveGroup, M::MatElem)
```

```@docs
representation_on_forms(G::LinearlyReductiveGroup, d::Int)
```

If `G` is a group returned by `linearly_reductive_group`, then

- `base_ring(G)` refers to $k$, and
- `defining_ideal(G)` to the ideal of the polynomials specifying $G$.

If `r` is a representation of `G`, then

- `group(r)` refers to $G$,
- `representation_matrix(r)` to the matrix specifying $\rho$, and
- `degree(r)` and `vector_space_dim(r)` to the dimension $n$ of $V$.

### Constructors for Invariant Rings

```@docs
invariant_ring(r::RepresentationLinearlyReductiveGroup)
```

## Basic Data Associated to Invariant Rings

If `RG` is the invariant ring $k[x_1, \dots, x_n]^G$ of a linearly reductive group $G$, then

- `group(RG)` refers to $G$,
- `coefficient_ring(RG)` to $k$, and
- `polynomial_ring(RG)` to $k[x_1, \dots, x_n]$.

## The Reynolds Operator

```@docs
reynolds_operator(R::RedGroupInvarRing, f::MPolyRingElem)
```

## The Null Cone

```@docs
null_cone_ideal(RG::RedGroupInvarRing)
```

## Fundamental Systems of Invariants

```@docs
fundamental_invariants(RG::RedGroupInvarRing, algorithm::Symbol = :default)
```

## Invariant Rings as Affine Algebras

```@docs
affine_algebra(RG::RedGroupInvarRing)
```
