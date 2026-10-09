# Linear Quotients

Linear quotients are orbit spaces of the action of a finite group $G$ on a finite-dimensional vector space $V$ over $\mathbb C$.
Formally, we define $V/G := \operatorname{Spec}\mathbb C[V]^G$.
Notice that the invariant ring $\mathbb C[V]^G$ is an affine algebra by a theorem of Hilbert-Noether.

## Constructor

Given a finite group $G\leq \operatorname{GL}_n(K)$, one can construct the
corresponding linear quotient $K^n/G$:
```@docs
linear_quotient(G::MatGroup)
```

!!! danger "Implicit choice of representation"
    Let $V = K^n$ be the regular representation of the matrix group $G$. In the current
    version, the object returned by `linear_quotient(G)` will work with the dual
    representation, that is, the linear quotient will be $V^\ast/G$. This might change
    in the future (notice that this code is still considered [experimental](https://docs.oscar-system.org/dev/Experimental/intro/))

!!! note "Root of unity"
    For many computations, we require that the base field `base_ring(G)` contains a
    primitive root of unity of order `exponent(G)`.
    If your chosen field is 'too small', you can easily change the base field with
    `map_entries(L, G)`, where `L` is the larger field.

## Class group

The divisor class group of a linear quotient $V/G$ is controlled by the pseudo-reflections
contained in the group $G$, see [Ben93](@cite).
```@docs
class_group(L::LinearQuotient)
```

## Singularities

One can study the types of the singularities of a linear quotient as follows.
```@docs
has_canonical_singularities(L::LinearQuotient)
```
```@docs
has_terminal_singularities(L::LinearQuotient)
```

## Cox rings of linear quotients

By a theorem of Arzhantsev and Gaifullin [AG10](@cite), the Cox ring of a linear
quotient $V/G$ is graded isomorphic to the invariant ring $K[V]^{[G,G]}$, where
$[G,G]$ is the derived subgroup of $G$.
```@docs
cox_ring(L::LinearQuotient)
```
