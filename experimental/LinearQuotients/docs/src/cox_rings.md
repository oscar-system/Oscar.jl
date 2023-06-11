```@meta
CurrentModule = Oscar
```

# Cox rings

## Cox rings of linear quotients

By a theorem of Arzhantsev and Gaifullin [AG10](@cite), the Cox ring of a linear
quotient $V/G$ is graded isomorphic to the invariant ring $K[V]^{[G,G]}$, where
$[G,G]$ is the derived subgroup of $G$.
This functionality is so far only available if the group does not contain any
reflections.
```@docs
cox_ring(L::LinearQuotient)
```

## Cox rings of $\mathbb Q$-factorial terminalizations

We provide an experimental algorithm to compute the Cox ring of a $\mathbb Q$-factorial
terminalization $X\to V/G$ of a linear quotient due to [Yam18](@cite).
```@docs
cox_ring_of_qq_factorial_terminalization(L::LinearQuotient)
```
