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
affine_hull(P::Polyhedron)
ambient_dim(P::Polyhedron)
dim(P::Polyhedron)
codim(P::Polyhedron)
isbounded(P::Polyhedron)
isfeasible(P::Polyhedron)
isfulldimensional(P::Polyhedron)
lineality_dim(P::Polyhedron)
lineality_space(P::Polyhedron)
recession_cone(P::Polyhedron)
relative_interior_point(P::Polyhedron)
```

## Combinatorial data

```@docs
nfacets(P::Polyhedron)
nvertices(P::Polyhedron)
f_vector(P::Polyhedron)
g_vector(P::Polyhedron)
h_vector(P::Polyhedron)
```

## Other

```@docs
boundary_lattice_points(P::Polyhedron)
contains(P::Polyhedron, v::AbstractVector)
ehrhart_polynomial(P::Polyhedron)
ehrhart_polynomial(R::FmpqPolyRing, P::Polyhedron)
h_star_polynomial(P::Polyhedron)
h_star_polynomial(R::FmpqPolyRing, P::Polyhedron)
interior_lattice_points(P::Polyhedron)
isnormal(P::Polyhedron)
issimple(P::Polyhedron)
issmooth(P::Polyhedron)
is_very_ample(P::Polyhedron)
lattice_points(P::Polyhedron)
lattice_volume(P::Polyhedron)
normalized_volume(P::Polyhedron)
polarize(P::Polyhedron)
project_full(P::Polyhedron)
print_constraints(A::AnyVecOrMat, b::AbstractVector; trivial::Bool = false)
print_constraints(P::Polyhedron; trivial::Bool = false)
solve_ineq(A::fmpz_mat, b::fmpz_mat)
solve_mixed(A::fmpz_mat, b::fmpz_mat, C::fmpz_mat, d::fmpz_mat)
solve_mixed(A::fmpz_mat, b::fmpz_mat, C::fmpz_mat)
solve_non_negative(A::fmpz_mat, b::fmpz_mat)
support_function(P::Polyhedron; convention::Symbol = :max)
volume(P::Polyhedron)
```


