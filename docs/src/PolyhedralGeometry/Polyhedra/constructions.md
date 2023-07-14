```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

# Constructions

The standard way to define a polyhedron is by either giving a
$V$-representation or an $H$-representation.  But polyhedra may also be
constructed through other means: by name, via operations on other polyhedra, or
from other objects in OSCAR.

## $H$- and $V$-representations

### Intersecting halfspaces: $H$-representation

```@docs
polyhedron(::Type{T}, A::AnyVecOrMat, b::AbstractVector) where T<:scalar_types
polyhedron(::Type{T}, I::Union{Nothing, AbstractCollection[AffineHalfspace]}, E::Union{Nothing, AbstractCollection[AffineHyperplane]} = nothing) where T<:scalar_types
```

The complete $H$-representation can be retrieved using [`facets`](@ref facets)
and [`affine_hull`](@ref affine_hull):
```jldoctest
julia> P = polyhedron(([-1 0; 1 0], [0,1]), ([0 1], [0]))
Polyhedron in ambient dimension 2

julia> facets(P)
2-element SubObjectIterator{AffineHalfspace{QQFieldElem}} over the Halfspaces of R^2 described by:
-x₁ ≦ 0
x₁ ≦ 1


julia> affine_hull(P)
1-element SubObjectIterator{AffineHyperplane{QQFieldElem}} over the Hyperplanes of R^2 described by:
x₂ = 0


julia> Q0 = polyhedron(facets(P))
Polyhedron in ambient dimension 2

julia> P == Q0
false

julia> Q1 = polyhedron(facets(P), affine_hull(P))
Polyhedron in ambient dimension 2

julia> P == Q1
true
```

### Computing convex hulls: $V$-representation

```@docs
convex_hull(::Type{T}, ::AnyVecOrMat; non_redundant::Bool=false) where T<:scalar_types
```

This is a standard triangle, defined via a (redundant) $V$-representation  and
its unique minimal $H$-representation:

```jldoctest
julia> T = convex_hull([ 0 0 ; 1 0 ; 0 1; 0 1/2 ])
Polyhedron in ambient dimension 2

julia> halfspace_matrix_pair(facets(T))
(A = [-1 0; 0 -1; 1 1], b = QQFieldElem[0, 0, 1])

```

The complete $V$-representation can be retrieved using [`minimal_faces`](@ref
minimal_faces), [`rays_modulo_lineality`](@ref rays_modulo_lineality) and [`lineality_space`](@ref lineality_space):

```jldoctest; filter = r"^polymake: +WARNING.*\n|^"
julia> P = convex_hull([0 0], [1 0], [0 1])
Polyhedron in ambient dimension 2

julia> Q0 = convex_hull(vertices(P))
Polyhedron in ambient dimension 2

julia> P == Q0
false

julia> mfP = minimal_faces(P)
(base_points = PointVector{QQFieldElem}[[0, 0]], lineality_basis = RayVector{QQFieldElem}[[0, 1]])

julia> rmlP = rays_modulo_lineality(P)
(rays_modulo_lineality = RayVector{QQFieldElem}[[1, 0]], lineality_basis = RayVector{QQFieldElem}[[0, 1]])

julia> Q1 = convex_hull(mfP.base_points, rmlP.rays_modulo_lineality)
Polyhedron in ambient dimension 2

julia> P == Q1
false

julia> Q0 == Q1
false

julia> Q2 = convex_hull(mfP.base_points, rmlP.rays_modulo_lineality, lineality_space(P))
Polyhedron in ambient dimension 2

julia> P == Q2
true
```

## Regular polytopes
A polytope is regular, in the strict sense, if it admits a flag-transtive group
of (linear) automorphisms. There are three infinite families of regular
polytopes which exist in each dimension: the (regular) simplices, cubes and
cross polytopes. In addition there are two exceptional regular 3-polytopes
(dodecahedron and icosahedron) plus three exceptional regular 4-polytopes
(24-cell, 120-cell and 600-cell).

The regular 3-polytopes are also known as the Platonic solids. Here we also
list the Archimedean, Catalan and Johnson solids, which form various
generalizations of the Platonic solids. However, here we implement "disjoint
families", i.e., the proper Archimedean solids exclude the Platonic solids;
similarly, the proper Johnson solids exclude the Archmidean solids.
```@docs
simplex
cross_polytope
cube
tetrahedron
dodecahedron
icosahedron
platonic_solid
archimedean_solid
johnson_solid
catalan_solid
regular_24_cell
regular_120_cell
regular_600_cell
```

## Other polytope constructions

```@docs
billera_lee_polytope
birkhoff_polytope
cyclic_polytope
del_pezzo_polytope
fano_simplex
fractional_cut_polytope
fractional_matching_polytope
gelfand_tsetlin_polytope
newton_polytope
orbit_polytope
rand_spherical_polytope
rand_subpolytope
```

## Operations on polyhedra
Polyhedra can be produced through operations on other polyhedra. For example,
they can be added using Minkowski addition or scaled; each of which results in
a new polyhedron.

```@docs
+(::Polyhedron{T}, ::Polyhedron{U}) where {T<:scalar_types, U<:scalar_types}
*(::Number, ::Polyhedron{T}) where T<:scalar_types
*(::Polyhedron{T}, ::Polyhedron{U})  where {T<:scalar_types, U<:scalar_types}
bipyramid
intersect(::Polyhedron...)
pyramid
```

The convex hull of two polytopes can be computed via `convex_hull`.
```@docs
convex_hull(::Polyhedron...)
```
