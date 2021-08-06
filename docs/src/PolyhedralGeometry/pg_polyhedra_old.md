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

Let $\mathbb{F}$ be an ordered field; the most prominent case here is where $\mathbb{F}=\mathbb{Q}$ are the rational numbers.

A set $P \subseteq \mathbb{F}^n$ is called a (convex) polyhedron if it can be written as the intersection of finitely many closed affine halfspaces in $\mathbb{F}^n$.
That is, there exists a matrix $A$ and a vector $b$ such that
$$P = P(A,b) = \{ x \in \mathbb{F}^n \mid Ax \leq b\}.$$
Writing $P$ as above is called an $H$-representation of $P$.

When a polyhedron $P \subset \mathbb{F}^n$ is bounded, it is called a polytope and the fundamental theorem of polytopes states that it may be written as the convex hull of finitely many points.
That is $$P = \textrm{conv}(p_1,\ldots,p_N), p_i \in \mathbb{F}^n.$$
Writing $P$ in this way is called a $V$-representation.
Polytopes are necessarily compact, i.e., they form convex bodies.

Each polytope has a unique $V$-representation which is minimal with respect to inclusion (or cardinality).
Conversely, a polyhedron which is full-dimensional, has a unique minimal $H$-representation.
If the polyhedron is not full-dimensional, then there is no canonical choice of an $H$-representation.


## Construction: H and V representations

Based on the definition of an $H$-representation, the constructor of `Polyhedron` can be used:

```@docs
Polyhedron(A::Union{Oscar.MatElem,AbstractMatrix}, b)
```

The following defines a polytope, as a `Polyhedron`, via a $V$-representation by calling the following function:

```@docs
convex_hull(::AnyVecOrMat; non_redundant::Bool=false)
```

### Computing convex hulls

This is a standard triangle, defined via a (redundant) $V$-representation  and its unique minimal $H$-representation:

```@repl oscar
T = convex_hull([ 0 0 ; 1 0 ; 0 1; 0 1/2 ])
facets_as_halfspace_matrix_pair(T)
```


## Other ways to construct polyhedra

Polyhedra may also be constructed by name, through operations on other polyhedra, or from other objects in Oscar.

### Constructing polyhedra by name
Some commonly used polyhedra are able to be constructed by name. For example

```@docs
simplex
```

```@docs
cube
```

```@docs
cross
```

```@docs
archimedean_solid
```

### Constructing polyhedra through operations
Polyhedra can be produced through operations on other polyhedra. For example, they can be added using Minkowski addition or scaled; each of which results in a new polyhedron.

```@docs
+(::Polyhedron, ::Polyhedron)
minkowski_sum
+(::AbstractVector, ::Polyhedron)
```

```@docs
*(::Int, ::Polyhedron)
```

You can also intersect polyhedra to obtain a new polyhedron.

```@docs
intersect(::Polyhedron, ::Polyhedron)
```


### Constructing polyhedra from other objects

There is a growing list of functionality for constructing polyhedra coming from other mathematical objects. We list a handful below.

```@docs
orbit_polytope
```

```@docs
newton_polytope
```





## `Polyhedron` and `polymake`'s `Polytope`

Many polyhedral computations are done through `polymake`.
This is visible in the structure `Polyhedron` via a pointer  `pm_polytope` to the corresponding `polymake` object.
This allows to apply all functionality of `polymake` by calling a suitable `Polymake.jl` function on this pointer.

To allow both `Oscar`'s and `polymake`'s functionality to be applicable to a polyhedron object, it can be converted back and forth:

```@docs
Polyhedron(::Polymake.BigObject)
pm_polytope(::Polyhedron)
```

There are several design differences between `polymake` and `Oscar`.

!!! warning
    Polyhedra in `polymake` and `Polymake.jl` use homogeneous coordinates. The polyhedra in `Oscar` use affine coordinates.
