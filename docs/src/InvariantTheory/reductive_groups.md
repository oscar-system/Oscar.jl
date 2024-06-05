```@meta
CurrentModule = Oscar
```

# Invariants of Linearly Reductive Groups

In this section, with notation as in the [introduction](@ref invariant_theory) to this chapter,
$G$ will be a *linearly algebraic group* over an algebraically closed
field $K$, $\rho: G \to \text{GL}(V)\cong \text{GL}_n(K)$ will
be a *rational* representation of $G$, and
$G$ will act on $K[V]\cong K[x_1, \dots, x_n]$ by linear
substitution: If $\rho(\pi) = (a_{i, j})$, then 

$(f \;\!   . \;\! \pi)  (x_1, \dots, x_n) = f\bigl(\sum_j a_{1, j}x_j, \dots, \sum_j a_{n, j}x_j\bigr).$

!!! note
    
    - The definition of linear reductivity guarantees the existence of a Reynolds operator $\mathcal R: K[V] \to K[V]$. 
    - By Hilbert's celebrated finiteness theorem, $K[V]^G$ is finitely generated as a $K$-algebra.
    - By a result of Hochster and Roberts, $K[V]^G$ is Cohen-Macaulay.

In cases where the Reynold's operator is explicitly known, generators of invariant rings of linearly reductive groups
can be found in two steps using Derksen's algorithm, see [Der99](@cite) :

- First, compute generators of Hilbert's null-cone ideal.
- Then, apply the Reynold's operator to these generators.

See also [DK15](@cite) and [DJ98](@cite).

## Creating Invariant Rings

### How Linearly Reductive Groups and Their Representations are Given

For the computation of invariant rings in the above setting, there is no need to deal with explicit elements of ``G`` or with its group structure.
The implementation of Derksen's algorithm in OSCAR can  handle situations where both $G$ and the representation $\rho$ are defined over an exact
subfield $k$ of $K$ which is supported by OSCAR: 

- ``G`` is  specified as an affine algebraic variety by polynomials with coefficients in $k$;
- ``\rho: G \to \text{GL}(V) \cong \text{GL}_n(K)`` is specified by an $n\times n$ matrix whose entries are polynomials in the same variables as those specifying $G$, with coefficients in $k$.

!!! note
    Proceeding as above is not a problem: Derksen's algorithms relies on Gröbner bases techniques and means to compute
    Reynolds operators. It does, thus, not change the initial ground field $k$. That is, all computations are performed over $k$
	and computations over any extension field of $k$ would lead to the same results.

In OSCAR, the basic set-up for a linearly reductive group in the context of Derksen's algorithm is provided by the 
function `linearly_reductive_group`. At current state, this only supports rational actions of  the special linear group
(in characteristic zero). For the action of this group by linear
substitution on, say, $n$-ary forms of degree $d$, an explicit Reynolds operator is
given by Cayley's Omega-process. We show this at work later in this section.


```@docs
linearly_reductive_group(sym::Symbol, m::Int, K::Field)
```

```@docs
linearly_reductive_group(sym::Symbol, m::Int, R::MPolyRing)
```

```@docs
representation_on_forms(G::LinearlyReductiveGroup, d::Int)
```

### Constructors for Invariant Rings

```@docs
invariant_ring(r::RepresentationLinearlyReductiveGroup)
```

```@docs
invariant_ring(R::MPolyDecRing, r::RepresentationLinearlyReductiveGroup)
```
 
 ### The Reynolds Operator

```@docs
reynolds_operator(R::RedGroupInvarRing, f::MPolyRingElem)
```

## Fundamental Systems of Invariants

```@docs
fundamental_invariants(RG::RedGroupInvarRing)
```





