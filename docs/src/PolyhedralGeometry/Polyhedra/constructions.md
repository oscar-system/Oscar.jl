# Constructions

The standard way to define a polyhedron is by either giving a
$V$-representation or an $H$-representation.  But polyhedra may also be
constructed through other means: by name, via operations on other polyhedra, or
from other objects in OSCAR.

## $H$- and $V$-representations

### Intersecting halfspaces: $H$-representation

```@docs
polyhedron(::Oscar.scalar_type_or_field, A::AnyVecOrMat, b::AbstractVector)
polyhedron(::Oscar.scalar_type_or_field, I::Union{Nothing, AbstractCollection[AffineHalfspace]}, E::Union{Nothing, AbstractCollection[AffineHyperplane]} = nothing)
```

The complete $H$-representation can be retrieved using [`facets`](@ref facets(as::Type{T}, P::Polyhedron{S}) where {S<:scalar_types,T<:Union{AffineHalfspace{S},AffineHyperplane{S},Pair{R,S} where R,Polyhedron{S}}})
and [`affine_hull`](@ref affine_hull(P::Polyhedron{T}) where {T<:scalar_types}):
```jldoctest
julia> P = polyhedron(([-1 0; 1 0], [0,1]), ([0 1], [0]))
Polyhedron in ambient dimension 2

julia> facets(P)
2-element SubObjectIterator{AffineHalfspace{QQFieldElem}} over the halfspaces of R^2 described by:
-x_1 <= 0
x_1 <= 1


julia> affine_hull(P)
1-element SubObjectIterator{AffineHyperplane{QQFieldElem}} over the hyperplanes of R^2 described by:
x_2 = 0


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
convex_hull(::Oscar.scalar_type_or_field, ::AnyVecOrMat; non_redundant::Bool=false)
```

This is a standard triangle, defined via a (redundant) $V$-representation  and
its unique minimal $H$-representation:

```jldoctest
julia> T = convex_hull([ 0 0 ; 1 0 ; 0 1; 0 1//2 ])
Polyhedron in ambient dimension 2

julia> halfspace_matrix_pair(facets(T))
(A = [-1 0; 0 -1; 1 1], b = QQFieldElem[0, 0, 1])

```

The complete $V$-representation can be retrieved using [`minimal_faces`](@ref
minimal_faces(P::Polyhedron{T}) where {T<:scalar_types}), [`rays_modulo_lineality`](@ref rays_modulo_lineality(P::Polyhedron{T}) where {T<:scalar_types}) and [`lineality_space`](@ref lineality_space(P::Polyhedron{T}) where {T<:scalar_types}):

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
A polytope is regular, in the strict sense, if it admits a flag-transitive group
of (linear) automorphisms. There are three infinite families of regular
polytopes which exist in each dimension: the (regular) simplices, cubes and
cross polytopes. In addition there are two exceptional regular 3-polytopes
(dodecahedron and icosahedron) plus three exceptional regular 4-polytopes
(24-cell, 120-cell and 600-cell).

The regular 3-polytopes are also known as the Platonic solids. Here we also
list the Archimedean, Catalan and Johnson solids, which form various
generalizations of the Platonic solids. However, here we implement "disjoint
families", i.e., the proper Archimedean solids exclude the Platonic solids;
similarly, the proper Johnson solids exclude the Archimedean solids.
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

Like some of the Johnson solids, the following four Archimedean and Catalan
solids are constructed using [serialized data](@ref serialization).
In order to properly document the respective sources, they also come as
separate functions.

```@docs
snub_cube
snub_dodecahedron
pentagonal_icositetrahedron
pentagonal_hexecontahedron
```

## Other polytope constructions

```@docs
SIM_body_polytope
associahedron
billera_lee_polytope
binary_markov_graph_polytope
birkhoff_polytope
cyclic_caratheodory_polytope
cyclic_polytope
del_pezzo_polytope
dwarfed_cube
dwarfed_product_polygons
explicit_zonotope
fano_simplex
fractional_cut_polytope
fractional_knapsack_polytope
fractional_matching_polytope
gelfand_tsetlin_polytope
goldfarb_cube
goldfarb_sit_cube
hypersimplex
hypertruncated_cube
k_cyclic_polytope
klee_minty_cube
lecture_hall_simplex
max_GC_rank_polytope
n_gon
newton_polytope
orbit_polytope
permutahedron
pile_polytope
pitman_stanley_polytope
perles_nonrational_8_polytope
pseudo_del_pezzo_polytope
rand01_polytope
rand_box_polytope
rand_cyclic_polytope
rand_metric
rand_metric_int
rand_normal_polytope
rand_spherical_polytope
rand_subpolytope
rss_associahedron
signed_permutahedron
stable_set_polytope
transportation_polytope
tutte_lifting
zonotope
zonotope_vertices_fukuda_matrix
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
vertex_figure
```

The convex hull of two polytopes can be computed via `convex_hull`.
```@docs
convex_hull(::Polyhedron...)
```
