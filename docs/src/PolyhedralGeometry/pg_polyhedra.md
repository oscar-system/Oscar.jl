```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["pg_polyhedra.md"]
```

# Polyhedra

## Introduction

A subset $P \subseteq \mathbb{R}^n$ is called a (metric) polyhedron if it can be written as the intersection of finitely many closed affine halfspaces in $\mathbb{R}^n$.
That is, there exists a matrix $A$ and a vector $b$ such that
$$P = P(A,b) = \{ x \in \mathbb{R}^n \mid Ax \leq b\}.$$
Writing $P$ as above is called the $H$-representation of $P$.

Polyhedra are convex bodies. When a polyhedron $P \subset \mathbb{R}^n$ is bounded, it is called a polytope and the fundamental theorem of polytopes states that it may be written as the convex hull of finitely many points in $\mathbb{R}^n$. That is $$P = \textrm{conv}(p_1,\ldots,p_N), p_i \in \mathbb{R}^n.$$ Writing $P$ in this way is called its $V$-representation.

Since many polyhedral computations are done through `polymake`, the structure `Polyhedron` contains a pointer  `pm_polytope` to the corresponding `polymake` object. Thus, all functionality of `Polymake.jl` is accessible by calling a `Polymake.jl` function on this pointer.

!!! warning
    Polyhedra in `polymake` and `Polymake.jl` use homogeneous coordinates. The polyhedra in `Oscar` use affine coordinates.

## Construction

Based on the definition of the $H$-representation, the constructor of `Polyhedron` can be used:

```@docs
Polyhedron(A::Union{Oscar.MatElem,AbstractMatrix}, b)
```

A `Polyhedron` given by its $V$-representation can also be created by calling the following function:

```@docs
convex_hull(::AnyVecOrMat)
```

Note that both examples provided result in equal polyhedra.

### Computing convex hulls

## `Polyhedron` and `polymake`'s `Polytope`

To allow both `Oscar`'s and `polymake`'s functionality to be applicable to a polyhedron object, it can be converted back and forth:

```@docs
Polyhedron(::Polymake.BigObject)
pm_polytope(::Polyhedron)
```
