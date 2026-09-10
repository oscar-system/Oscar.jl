# Local Cohomology
Let $M$ be a finitely generated $\mathbb{Z}^d$-graded module over a monoid algebra $k[Q]$. 
Further, let $I\subseteq k[Q]$ be an ideal. 
The $i$-th local cohomology module of $M$, denoted $H^i_I(M)$, is obtained as follows:

Let

$J^\bullet \colon 0 \to M \xrightarrow{\epsilon} J^0 \xrightarrow{d^0} J^1 \xrightarrow{d^1}\cdots \xrightarrow{d^{i-1}} J^i \xrightarrow{d^i} \cdots$

be an injective resolution of $M$. Applying the left exact functor $\Gamma_I$, which maps a $\mathbb{Z}^d$-graded module $N$ to the submodule

$\Gamma_I(N) = \{n \in N \mid \exists m \in \mathbb{N} \colon n\cdot I^m = 0\},$

to $J^\bullet$ we obtain the complex

$\Gamma_I(J^\bullet) \colon 0 \to \Gamma_I(J^0) \xrightarrow{d^0} \Gamma_I(J^1) \xrightarrow{d^1}\cdots \xrightarrow{d^{i-1}} \Gamma_I(J^i) \xrightarrow{d^i} \cdots.$

The *$i$-th local cohomology module of $M$ supported on $I$* is the $i$-th cohomology module of $\Gamma_I(J^\bullet)$.

!!! note
    We require that the monoid algebra $k[Q]$ is normal. 

## Cohomological degree zero
The zeroth local cohomology module of $M$ supported on $I$ is

$H^0_I(M) = \Gamma_I(M) = \{m\in M \mid \exists n \in \mathbb{N} \colon m\cdot I^n = 0\}.$

The function `zeroth_local_cohomology` returns this submodule.
The function is split off from the other degrees because it returns its result 
as a proper module in OSCAR while higher local cohomology modules are
represented using sector partitions.

```@docs
zeroth_local_cohomology(M::SubquoModule{<:MonoidAlgebraElem}, I::MonoidAlgebraIdeal)
```

## Sector partitions of local cohomology modules
The local cohomology modules $H^i_I(M)$ are in general not finitely generated 
for $i>0$. However, sector partitions are a finite data structure for them.

A *sector partition* $\mathcal{S}$ of $H^i_I(M)$ consists of

- a finite partition $\mathbb{Z}^d = \sqcup_{S\in \mathcal{S}} S$ into *sectors*
- finite dimensional $k$-vector spaces $H_S$ for each sector $S\in \mathcal{S}$, and,
- maps between these vector spaces.

Now given $\alpha \in \mathbb{Z}^d$,

$H^i_I(M)_\alpha \cong k^{\dim(H_S)} \text{ for } \alpha \in S.$

For more details on sector partitions see, e.g., Chapter 13 of [MS05](@cite).

The function [`local_cohomology`](@ref) computes a sector partition 
of $H^i_I(M)$. For performance, multiple local cohomology modules $H^1_I(M),\dots,H^i_I(M)$ 
should be computed at once using the function `local_cohomology_all`.

```@docs
    local_cohomology(M::SubquoModule{T}, I::MonoidAlgebraIdeal, i::Integer) where {T<:MonoidAlgebraElem}
    local_cohomology_all(M::SubquoModule{T}, I::MonoidAlgebraIdeal, i::Integer) where {T<:MonoidAlgebraElem}
```

### Data associated to sector partitions
Let `H = local_cohomology(M, I, i)` be a sector partition of the local cohomology module $H^i_I(M)$.
Then `sectors(H)` returns the finite partition of $\mathbb{Z}^d$ into sectors. Each sector `S` consists of

- the finite dimensional $k$-vector space $H_S$ = `S.H`, and
- the sector as a polyhedron `S.sector`.

### Tests on local cohomology modules
To test vanishing of local cohomology independent of the internal representation use `is_zero`:
```@docs
is_zero(S::Oscar.InjectiveResolutions.SectorPartitionLC)
```
