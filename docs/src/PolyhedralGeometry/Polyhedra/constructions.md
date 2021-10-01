```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["constructions.md"]
```

# Constructions

The standard way to define a polyhedron is by either giving a
$V$-representation or an $H$-representation.  But polyhedra may also be
constructed through other means: by name, via operations on other polyhedra, or
from other objects in Oscar.

## $H$- and $V$-representations

### Intersecting halfspaces: $H$-representation


```@docs
Polyhedron(A::Union{Oscar.MatElem,AbstractMatrix}, b)
```

### Computing convex hulls: $V$-representation


```@docs
convex_hull(::AnyVecOrMat; non_redundant::Bool=false)
```

This is a standard triangle, defined via a (redundant) $V$-representation  and
its unique minimal $H$-representation:

```@repl oscar
T = convex_hull([ 0 0 ; 1 0 ; 0 1; 0 1/2 ])
halfspace_matrix_pair(facets(T))
```




## Named polyhedra

```@docs
simplex
cube
cross
archimedean_solid
```

## Operations on polyhedra
Polyhedra can be produced through operations on other polyhedra. For example,
they can be added using Minkowski addition or scaled; each of which results in
a new polyhedron.

```@docs
+(::Polyhedron, ::Polyhedron)
*(::Int, ::Polyhedron)
intersect(::Polyhedron, ::Polyhedron)
*(::Polyhedron, ::Polyhedron)
```

The convex hull of two polytopes can be computed via `convex_hull`.
```@docs
convex_hull(::Polyhedron,::Polyhedron)
```



## Polyhedra from other mathematical objects


```@docs
orbit_polytope
newton_polytope
```
