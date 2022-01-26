```@meta
CurrentModule = Oscar
DocTestSetup = quote
   using Oscar
end
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

The complete $H$-representation can be retrieved using [`facets`](@ref facets)
and [`affine_hull`](@ref affine_hull):
```jldoctest
julia> P = Polyhedron(([-1 0; 1 0], [0,1]), ([0 1], [0]))
A polyhedron in ambient dimension 2

julia> facets(P)
2-element SubObjectIterator{AffineHalfspace}:
 The Halfspace of R^2 described by
1: -x₁ ≦ 0

 The Halfspace of R^2 described by
1: x₁ ≦ 1


julia> affine_hull(P)
1-element SubObjectIterator{AffineHyperplane}:
 The Hyperplane of R^2 described by
1: x₂ = 0


julia> Q0 = Polyhedron(facets(P))
A polyhedron in ambient dimension 2

julia> P == Q0
false

julia> Q1 = Polyhedron(facets(P), affine_hull(P))
A polyhedron in ambient dimension 2

julia> P == Q1
true
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

The complete $V$-representation can be retrieved using [`vertices`](@ref
vertices), [`rays`](@ref rays) and [`lineality_space`](@ref lineality_space):

```jldoctest; filter = r"^polymake: +WARNING.*\n|^"
julia> P = convex_hull([0 0], [1 0], [0 1])
A polyhedron in ambient dimension 2

julia> Q0 = convex_hull(vertices(P))
A polyhedron in ambient dimension 2

julia> P == Q0
false

julia> Q1 = convex_hull(vertices(P), rays(P))
A polyhedron in ambient dimension 2

julia> P == Q1
false

julia> Q0 == Q1
false

julia> Q2 = convex_hull(vertices(P), rays(P), lineality_space(P))
A polyhedron in ambient dimension 2

julia> P == Q2
true
```

## Named polyhedra

```@docs
archimedean_solid
birkhoff
cross
cube
gelfand_tsetlin
simplex
```

## Operations on polyhedra
Polyhedra can be produced through operations on other polyhedra. For example,
they can be added using Minkowski addition or scaled; each of which results in
a new polyhedron.

```@docs
+(::Polyhedron, ::Polyhedron)
*(::Int, ::Polyhedron)
*(::Polyhedron, ::Polyhedron)
bipyramid
intersect(::Polyhedron, ::Polyhedron)
pyramid
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
