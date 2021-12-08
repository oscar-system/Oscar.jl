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

```@docs
affine_hull(P::Polyhedron)
ambient_dim(P::Polyhedron)
boundary_lattice_points(P::Polyhedron)
codim(P::Polyhedron)
contains(P::Polyhedron, v::AbstractVector)
dim(P::Polyhedron)
ehrhart_polynomial(P::Polyhedron)
ehrhart_polynomial(R::FmpqMPolyRing, P::Polyhedron)
facets(P::Polyhedron)
f_vector(P::Polyhedron)
interior_lattice_points(P::Polyhedron)
isbounded(P::Polyhedron)
isfeasible(P::Polyhedron)
isfulldimensional(P::Polyhedron)
isnormal(P::Polyhedron)
issimple(P::Polyhedron)
issmooth(P::Polyhedron)
isvery_ample(P::Polyhedron)
lattice_points(P::Polyhedron)
lattice_volume(P::Polyhedron)
lineality_dim(P::Polyhedron)
lineality_space(P::Polyhedron)
nfacets(P::Polyhedron)
normalized_volume(P::Polyhedron)
nvertices(P::Polyhedron)
polarize(P::Polyhedron)
project_full(P::Polyhedron)
print_constraints(A::AnyVecOrMat, b::AbstractVector; trivial::Bool = false)
print_constraints(P::Polyhedron; trivial::Bool = false)
rays(P::Polyhedron)
recession_cone(P::Polyhedron)
relative_interior_point(P::Polyhedron)
solve_ineq(A::fmpz_mat, b::fmpz_mat)
solve_mixed(A::fmpz_mat, b::fmpz_mat, C::fmpz_mat, d::fmpz_mat)
solve_mixed(A::fmpz_mat, b::fmpz_mat, C::fmpz_mat)
solve_non_negative(A::fmpz_mat, b::fmpz_mat)
support_function(P::Polyhedron; convention::Symbol = :max)
vertices(P::Polyhedron)
volume(P::Polyhedron)
```


