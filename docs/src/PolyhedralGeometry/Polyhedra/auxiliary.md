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
ambient_dim(P::Polyhedron)
codim(P::Polyhedron)
contains(P::Polyhedron, v::AbstractVector)
dim(P::Polyhedron)
facets(P::Polyhedron)
f_vector(P::Polyhedron)
isbounded(P::Polyhedron)
isfeasible(P::Polyhedron)
isfulldimensional(P::Polyhedron)
isnormal(P::Polyhedron)
issmooth(P::Polyhedron)
lattice_points(P::Polyhedron)
lineality_dim(P::Polyhedron)
lineality_space(P::Polyhedron)
nfacets(P::Polyhedron)
normalized_volume(P::Polyhedron)
nvertices(P::Polyhedron)
polarize(P::Polyhedron)
print_constraints(A::AnyVecOrMat, b::AbstractVector; trivial::Bool = false)
print_constraints(P::Polyhedron; trivial::Bool = false)
rays(P::Polyhedron)
recession_cone(P::Polyhedron)
support_function(P::Polyhedron; convention::Symbol = :max)
volume(P::Polyhedron)
```


