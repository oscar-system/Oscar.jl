```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Local Cohomology
Let $M$ be a finitely generated $\mathbb{Z}^d$-graded module over a monoid algebra $k[Q]$. Further, let $J\subseteq k[Q]$ be an ideal. The $i$-th local cohomology module of $M$, denoted $H^i_J(M)$, is obtained as follows:

Let

$I^\bullet \colon 0 \to M \xrightarrow{\epsilon} I^0 \xrightarrow{d^0} I^1 \xrightarrow{d^1}\cdots \xrightarrow{d^{i-1}} I^i \xrightarrow{d^i} \cdots$

be an injective resolution of $M$. Applying the left exact functor $\Gamma_J$, that maps a $\mathbb{Z}^d$-graded module $N$ to the submodule

$\Gamma_I(N) = \{n \in \mid \exists n \in \mathbb{N} \colon n\cdot J^n = 0\},$

to $I^\bullet$ we obtain the complex

$\Gamma_J(I^\bullet) \colon 0 \to \Gamma_J(I^0) \xrightarrow{d^0} \Gamma_J(I^1) \xrightarrow{d^1}\cdots \xrightarrow{d^{i-1}} \Gamma_J(I^i) \xrightarrow{d^i} \cdots.$

The *$i$-th local cohomology module of $M$ supported on $J$* is the $i$-th cohomology module of $\Gamma_J(I^\bullet)$.

!!! note
    We require that the monoid algebra $k[Q]$ is normal. 

## $H^0_J(M)$

## Computing Local Cohomology


### Data asssociated to local cohomology


```@docs
zeroth_local_cohomology

local_cohomology

local_cohomology_all
```