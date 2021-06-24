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

A set $P \subseteq \mathbb{R}^n$ is called a (convex) polyhedron if it can be written as the intersection of finitely many closed affine halfspaces in $\mathbb{R}^n$.
That is, there exists a matrix $A$ and a vector $b$ such that
$$P = P(A,b) = \{ x \in \mathbb{R}^n \mid Ax \leq b\}.$$
Writing $P$ as above is called an $H$-representation of $P$.

When a polyhedron $P \subset \mathbb{R}^n$ is bounded, it is called a polytope and the fundamental theorem of polytopes states that it may be written as the convex hull of finitely many points.
That is $$P = \textrm{conv}(p_1,\ldots,p_N), p_i \in \mathbb{R}^n.$$
Writing $P$ in this way is called a $V$-representation.
Polytopes are necessarily compact, i.e., they form convex bodies.

Each polytope has a unique $V$-representation which is minimal with respect to inclusion (or cardinality).
Conversely, a polyhedron which is full-dimensional, has a unique minimal $H$-representation.
If the polyhedron is not full-dimensional, then there is no canonical choice of an $H$-representation.

Since many polyhedral computations are done through `polymake`, the structure `Polyhedron` contains a pointer  `pm_polytope` to the corresponding `polymake` object.
Thus, all functionality of `Polymake.jl` is accessible by calling a suitable `Polymake.jl` function on this pointer.

!!! warning
    Polyhedra in `polymake` and `Polymake.jl` use homogeneous coordinates. The polyhedra in `Oscar` use affine coordinates.

## Construction

Based on the definition of an $H$-representation, the constructor of `Polyhedron` can be used:

```@docs
Polyhedron(A::Union{Oscar.MatElem,AbstractMatrix}, b)
```

The following defines a polytope, as a `Polyhedron`, via a $V$-representation by calling the following function:

```@docs
convex_hull(::AnyVecOrMat)
```

### Computing convex hulls

This is a standard triangle, defined via a (redundant) $V$-representation  and its unique minimal $H$-representation:

```@repl oscar
T = convex_hull([ 0 0 ; 1 0 ; 0 1; 0 1/2 ])
facets_as_halfspace_matrix_pair(T)
```

## `Polyhedron` and `polymake`'s `Polytope`

To allow both `Oscar`'s and `polymake`'s functionality to be applicable to a polyhedron object, it can be converted back and forth:

```@docs
Polyhedron(::Polymake.BigObject)
pm_polytope(::Polyhedron)
```
