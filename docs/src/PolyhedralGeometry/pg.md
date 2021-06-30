```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["pg.md"]
```

# Introduction
```@docs
Halfspaces
faces
vertices
vertices_as_point_matrix
nrays
nvertices
rays
nfacets
facets
facets_as_halfspace_matrix_pair
volume
normalized_volume
dim
lattice_points
ambient_dim
codim
lineality_space
recession_cone
isfeasible
issmooth
isnormal
isbounded
isfulldimensional
f_vector
support_function
orbit_polytope
cube
newton_polytope
intersect
minkowski_sum
+(::Polyhedron, ::Polyhedron)
*(::Polyhedron, ::Int)
*(::Int, ::Polyhedron)
+(::Polyhedron, ::AbstractVector)
+(::AbstractVector, ::Polyhedron)
cross
archimedean_solid
```

The polyhedral geometry part of OSCAR provides functionality for handling
- convex polytopes, unbounded polyhedra and cones
- linear programs

General textbooks offering details on theory and algorithms include:
- [JT13](@cite)
- [Zie95](@cite)
