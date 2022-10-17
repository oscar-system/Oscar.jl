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
Polyhedron{T}(A::AnyVecOrMat, b::AbstractVector) where T<:scalar_types
Polyhedron{T}(I::Union{Nothing, AbstractCollection[AffineHalfspace]}, E::Union{Nothing, AbstractCollection[AffineHyperplane]} = nothing) where T<:scalar_types
```

The complete $H$-representation can be retrieved using [`facets`](@ref facets)
and [`affine_hull`](@ref affine_hull):
```jldoctest
julia> P = Polyhedron(([-1 0; 1 0], [0,1]), ([0 1], [0]))
A polyhedron in ambient dimension 2

julia> facets(P)
2-element SubObjectIterator{AffineHalfspace{fmpq}} over the Halfspaces of R^2 described by:
-x₁ ≦ 0
x₁ ≦ 1


julia> affine_hull(P)
1-element SubObjectIterator{AffineHyperplane{fmpq}} over the Hyperplanes of R^2 described by:
x₂ = 0


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
convex_hull(::Type{T}, ::AnyVecOrMat; non_redundant::Bool=false) where T<:scalar_types
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
catalan_solid
cross
cube
cyclic_polytope
del_pezzo_polytope
fano_simplex
fractional_cut_polytope
fractional_matching_polytope
gelfand_tsetlin
simplex
```

## Operations on polyhedra
Polyhedra can be produced through operations on other polyhedra. For example,
they can be added using Minkowski addition or scaled; each of which results in
a new polyhedron.

```@docs
+(::Polyhedron{T}, ::Polyhedron{T}) where T<:scalar_types
*(::Int, ::Polyhedron{T}) where T<:scalar_types
*(::Polyhedron{T}, ::Polyhedron{T})  where T<:scalar_types
bipyramid
intersect(::Polyhedron{T}, ::Polyhedron{T}) where T<:scalar_types
pyramid
```

The convex hull of two polytopes can be computed via `convex_hull`.
```@docs
convex_hull(::Polyhedron{T},::Polyhedron{T}) where T<:scalar_types
```

## Polyhedra from other mathematical objects


```@docs
orbit_polytope
newton_polytope
```
