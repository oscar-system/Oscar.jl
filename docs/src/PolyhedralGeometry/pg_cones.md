```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["pg_cones.md"]
```

# Cones

## Introduction

A set $C \subseteq \mathbb{R}^n$ is called a (polyhedral) cone if it can be written as the intersection of finitely many closed (non-affine) halfspaces in $\mathbb{R}^n$. This especially means that the origin is contained in each of these halfspaces.
Looking at that statement from an algebraic point of view, there exists a matrix $A$ such that
$$C = C(A) = \{ x \in \mathbb{R}^n \mid Ax \leq 0\}.$$
Writing $C$ as above is called an $H$-representation of $C$.

A cone $C \subset \mathbb{R}^n$ may also be written as the positive hull of finitely many points.
That is $$C = \textrm{pos}(p_1,\ldots,p_N) = \{ \sum_{i = 1}^N a_i p_i\ |\ a_i \in \mathbb{R}_{\geq 0} \}, p_i \in \mathbb{R}^n.$$
Writing $C$ in this way is called a $V$-representation.
Cones are necessarily convex.

Each cone has a unique $V$-representation which is minimal with respect to inclusion (or cardinality).
Conversely, a cone which is full-dimensional, has a unique minimal $H$-representation.
If the cone is not full-dimensional, then there is no canonical choice of an $H$-representation.

## Constructions

The following function defines a `Cone` via a $V$-representation:

```@docs
convex_hull(::AnyVecOrMat)
```
