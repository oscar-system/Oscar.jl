```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["auxiliary.md"]
```

# Auxiliary functions

## Geometric data

```@docs
facets(P::Polyhedron)
vertices(P::Polyhedron)
rays(P::Polyhedron)
affine_hull(P::Polyhedron{T}) where T<:scalar_types
ambient_dim(P::Polyhedron)
dim(P::Polyhedron)
codim(P::Polyhedron)
is_bounded(P::Polyhedron)
is_feasible(P::Polyhedron)
is_fulldimensional(P::Polyhedron)
lineality_dim(P::Polyhedron)
lineality_space(P::Polyhedron{T}) where T<:scalar_types
recession_cone(P::Polyhedron{T}) where T<:scalar_types
relative_interior_point(P::Polyhedron{T}) where T<:scalar_types
```

## Combinatorial data

```@docs
nfacets(P::Polyhedron)
nvertices(P::Polyhedron)
f_vector(P::Polyhedron)
g_vector(P::Polyhedron)
h_vector(P::Polyhedron)
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
boundary_lattice_points(P::Polyhedron{fmpq})
contains(P::Polyhedron, v::AbstractVector)
ehrhart_polynomial(P::Polyhedron{fmpq})
ehrhart_polynomial(R::FmpqPolyRing, P::Polyhedron{fmpq})
h_star_polynomial(P::Polyhedron{fmpq})
h_star_polynomial(R::FmpqPolyRing, P::Polyhedron{fmpq})
interior_lattice_points(P::Polyhedron{fmpq})
is_normal(P::Polyhedron{fmpq})
is_simple(P::Polyhedron)
is_smooth(P::Polyhedron{fmpq})
is_very_ample(P::Polyhedron{fmpq})
lattice_points(P::Polyhedron{fmpq})
lattice_volume(P::Polyhedron{fmpq})
normalized_volume(P::Polyhedron{T}) where T<:scalar_types
polarize(P::Polyhedron{T}) where T<:scalar_types
project_full(P::Polyhedron{T}) where T<:scalar_types
print_constraints(A::AnyVecOrMat, b::AbstractVector; trivial::Bool = false)
print_constraints(P::Polyhedron; trivial::Bool = false)
regular_triangulations
regular_triangulation
secondary_polytope
solve_ineq(A::fmpz_mat, b::fmpz_mat)
solve_mixed(A::fmpz_mat, b::fmpz_mat, C::fmpz_mat, d::fmpz_mat)
solve_mixed(A::fmpz_mat, b::fmpz_mat, C::fmpz_mat)
solve_non_negative(A::fmpz_mat, b::fmpz_mat)
support_function(P::Polyhedron{T}; convention::Symbol = :max) where T<:scalar_types
volume(P::Polyhedron{T}) where T<:scalar_types
```


