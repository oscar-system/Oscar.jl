```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["pg_polyhedra_constructions.md"]
```

# Constructions

The standard way to define a polyhedron is by either giving a $V$-representation or an $H$-representation.
But polyhedra may also be constructed through other means: by name, via operations on other polyhedra, or from other objects in Oscar.

## H and V representations

### Intersecting halfspaces: H-representation

Based on the definition of an $H$-representation, the constructor of `Polyhedron` can be used:

```@docs
Polyhedron(A::Union{Oscar.MatElem,AbstractMatrix}, b)
```

The following defines a polytope, as a `Polyhedron`, via a $V$-representation by calling the following function:

```@docs
convex_hull(::AnyVecOrMat; non_redundant::Bool=false)
```

### Computing convex hulls: V-representation

This is a standard triangle, defined via a (redundant) $V$-representation  and its unique minimal $H$-representation:

```@repl oscar
T = convex_hull([ 0 0 ; 1 0 ; 0 1; 0 1/2 ])
facets_as_halfspace_matrix_pair(T)
```




## Named polyhedra
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

## Operations on polyhedra
Polyhedra can be produced through operations on other polyhedra. For example, they can be added using Minkowski addition or scaled; each of which results in a new polyhedron.

```@docs
+(::Polyhedron, ::Polyhedron)
```

```@docs
*(::Int, ::Polyhedron)
```

You can also intersect polyhedra to obtain a new polyhedron, or take Cartesian products.

```@docs
intersect(::Polyhedron, ::Polyhedron)
```

```@docs
*(::Polyhedron, ::Polyhedron)
```

The convex hull of two polytopes is supported by the function `convex_hull`.
```@docs
convex_hull(::Polyhedron,::Polyhedron)
```



## Polyhedra from other objects

There is a growing list of functionality for constructing polyhedra coming from other mathematical objects. We list a handful below.

```@docs
orbit_polytope
```

```@docs
newton_polytope
```
