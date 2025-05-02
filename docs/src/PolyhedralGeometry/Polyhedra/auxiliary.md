```@meta
CurrentModule = Oscar
```

# Auxiliary functions

## Geometric data

```@docs
facets(as::Type{T}, P::Polyhedron{S}) where {S<:scalar_types,T<:Union{AffineHalfspace{S},AffineHyperplane{S},Pair{R,S} where R,Polyhedron{S}}}
vertices(as::Type{PointVector{T}}, P::Polyhedron{T}) where {T<:scalar_types}
rays(as::Type{RayVector{T}}, P::Polyhedron{T}) where {T<:scalar_types}
rays_modulo_lineality(P::Polyhedron{T}) where T<:scalar_types
minimal_faces(P::Polyhedron{T}) where T<:scalar_types
affine_hull(P::Polyhedron{T}) where T<:scalar_types
ambient_dim(P::Polyhedron)
dim(P::Polyhedron)
codim(P::Polyhedron)
is_bounded(P::Polyhedron)
is_feasible(P::Polyhedron)
is_fulldimensional(P::Polyhedron)
is_lattice_polytope(P::Polyhedron{QQFieldElem})
lineality_dim(P::Polyhedron)
lineality_space(P::Polyhedron{T}) where T<:scalar_types
recession_cone(P::Polyhedron{T}) where T<:scalar_types
relative_interior_point(P::Polyhedron{T}) where T<:scalar_types
```

## Combinatorial data

```@docs
n_facets(P::Polyhedron)
n_vertices(P::Polyhedron)
f_vector(P::Polyhedron)
facet_sizes(P::Polyhedron)
g_vector(P::Polyhedron)
h_vector(P::Polyhedron)
vertex_sizes(P::Polyhedron)
```

## Groups
```@docs
linear_symmetries(P::Polyhedron)
combinatorial_symmetries(P::Polyhedron)
automorphism_group(P::Polyhedron)
automorphism_group_generators(P::Polyhedron)
automorphism_group(IM::IncidenceMatrix)
automorphism_group_generators(IM::IncidenceMatrix)
```

## Other

```@docs
all_triangulations
boundary_lattice_points(P::Polyhedron{QQFieldElem})
Base.in(v::AbstractVector, P::Polyhedron)
Base.issubset(P::Polyhedron{T}, Q::Polyhedron{T}) where T<:scalar_types
demazure_character(lambda::AbstractVector, sigma::PermGroupElem)
ehrhart_polynomial(P::Polyhedron{QQFieldElem})
ehrhart_polynomial(R::QQPolyRing, P::Polyhedron{QQFieldElem})
h_star_polynomial(P::Polyhedron{QQFieldElem})
h_star_polynomial(R::QQPolyRing, P::Polyhedron{QQFieldElem})
interior_lattice_points(P::Polyhedron{QQFieldElem})
is_normal(P::Polyhedron{QQFieldElem})
is_simple(P::Polyhedron)
is_smooth(P::Polyhedron{QQFieldElem})
is_very_ample(P::Polyhedron{QQFieldElem})
is_archimedean_solid(P::Polyhedron)
is_johnson_solid(P::Polyhedron)
is_platonic_solid(P::Polyhedron)
lattice_points(P::Polyhedron{QQFieldElem})
lattice_volume(P::Polyhedron{QQFieldElem})
normalized_volume(P::Polyhedron)
polarize(P::Polyhedron{T}) where T<:scalar_types
project_full(P::Polyhedron{T}) where T<:scalar_types
print_constraints(io::IO, A::AnyVecOrMat, b::AbstractVector; trivial::Bool = false, numbered::Bool = false, cmp = :lte)
print_constraints(io::IO, P::Polyhedron; trivial::Bool = false, numbered::Bool = false)
regular_triangulations
regular_triangulation
secondary_polytope
solve_ineq(as::Type{T}, A::ZZMatrix, b::ZZMatrix) where {T}
solve_mixed(as::Type{T}, A::ZZMatrix, b::ZZMatrix, C::ZZMatrix, d::ZZMatrix) where {T}
solve_mixed(as::Type{T}, A::ZZMatrix, b::ZZMatrix, C::ZZMatrix) where {T}
solve_non_negative(as::Type{T}, A::ZZMatrix, b::ZZMatrix) where {T}
support_function(P::Polyhedron{T}; convention::Symbol = :max) where T<:scalar_types
volume(P::Polyhedron{T}) where T<:scalar_types
```


