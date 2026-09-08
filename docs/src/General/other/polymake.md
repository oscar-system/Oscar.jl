# Notes for polymake users

OSCAR uses polymake for polyhedral geometry and parts of combinatorics,
via the Julia package
[Polymake.jl](https://github.com/oscar-system/Polymake.jl),
which is available as `Polymake`.
The first section describes how to reach polymake functionality that has
no OSCAR counterpart yet, the remaining sections describe differences
between polymake and OSCAR, and list common polymake commands together
with their OSCAR counterparts.

!!! note "Help wanted"
    This page is a start. Please tell us what is missing, see
    [Notes for users of other computer algebra systems](@ref).

## Calling polymake from OSCAR

Polyhedra, cones, fans, matroids, simplicial complexes, and graphs in
OSCAR wrap polymake objects, and `Oscar.pm_object` returns the wrapped
object.
Properties of polymake objects are accessed with the dot syntax,
and the functions of a polymake application `app` are available as
`Polymake.app.<name>`.
The results are polymake objects; `polyhedron` turns a polymake polytope
back into an OSCAR polyhedron, and `matrix` converts polymake matrices.
```jldoctest
julia> P = cube(3);

julia> pm = Oscar.pm_object(P);

julia> pm.F_VECTOR
pm::Vector<pm::Integer>
8 12 6

julia> matrix(QQ, pm.VERTICES)
[1   -1   -1   -1]
[1    1   -1   -1]
[1   -1    1   -1]
[1    1    1   -1]
[1   -1   -1    1]
[1    1   -1    1]
[1   -1    1    1]
[1    1    1    1]

julia> polyhedron(Polymake.polytope.hypersimplex(2, 4))
Polytope in ambient dimension 4
```
Note that `Oscar.pm_object` is not part of OSCAR's official interface
and may change.
Prefer the OSCAR functions whenever they exist
(`f_vector(P)`, `vertices(P)`, and `hypersimplex(2, 4)` in the example
above), they take care of the differences described below, and please
[open an issue](https://github.com/oscar-system/Oscar.jl/issues)
if you need polymake functionality that OSCAR does not provide.
More details can be found in the
[Polymake.jl manual](https://oscar-system.github.io/Polymake.jl/stable/).

## Differences in syntax

- Polymake's language is Perl, so variables start with `$`,
  every statement ends with `;`, and values are shown with `print`.
  In Julia, `P = cube(3)` assigns, and entering `P` or `f_vector(P)`
  shows the value, see [Semicolons and output](@ref).

- Properties become functions:
  `$p->VERTICES` is `vertices(P)`, `$p->N_VERTICES` is `n_vertices(P)`,
  `$p->SIMPLE` is `is_simple(P)`.
  Objects are constructed by functions, too:
  `new Polytope(POINTS=>[[1,0,0],[1,1,0],[1,0,1]])` is
  `convex_hull([0 0; 1 0; 0 1])`,
  and `new Cone(INPUT_RAYS=>[[1,0],[1,1]])` is `positive_hull([1 0; 1 1])`.

- `user_method`s cannot be accessed via Julia's dot syntax, i.e. something like

  ```julia
  c = Polymake.polytope.cube(3)
  c.AMBIENT_DIM
  ```

  will not work. Instead `user_method`s are attached as Julia functions in
  their respective application. They are always written in lowercase. In the
  example the following works:

  ```julia
  c = Polymake.polytope.cube(3)
  Polymake.polytope.ambient_dim(c)
  ```

- There are no applications to switch between.
  `application "fan";` is not needed, `normal_fan(P)` is available
  directly, and so are matroids, graphs, and simplicial complexes.

- `help "cube";` is `?cube`, and `$p->VISUAL;` is `visualize(P)`.

- `save($p, "file.poly");` and `load("file.poly")` are
  `save("file.mrdi", P)` and `load("file.mrdi")`,
  see [Serialization](@ref).

## Differences in semantics

- OSCAR (and Julia) is `1`-based, meaning that it counts from `1`, rather than
  from `0` like polymake. For most properties we have taken care of the
  translation but be aware that it might pop up at some point and generate
  confusion.

  For convenience, `Polymake.jl` provides `Polymake.to_one_based_indexing` and
  `Polymake.to_zero_based_indexing`.

- Polyhedra and polyhedral complexes in OSCAR are represented inhomogeneously,
  i.e. without the leading `1` for vertices or `0` for rays. Hence constructors
  take points, rays, and lineality generators separately.
  Similarly, a row `[b, a1, a2]` of polymake's `INEQUALITIES` encodes
  ``b + a_1 x_1 + a_2 x_2 \geq 0``,
  whereas OSCAR's `polyhedron(A, b)` encodes ``Ax \leq b``.
  The homogeneous coordinates are still visible when accessing the
  polymake object directly, as in the example above.

- Properties are computed on demand in both systems.
  A polymake object stores every property that has been computed;
  the OSCAR functions store their results in the wrapped polymake
  object in the same way, so calling `volume(P)` twice computes the
  volume only once.

- Many OSCAR functions return iterators over OSCAR objects instead of
  matrices, for example `vertices(P)` yields the vertices as vectors,
  `facets(P)` yields halfspaces, and `faces(P, 1)` yields the edges as
  polyhedra.
  Use `point_matrix(vertices(P))` for the matrix of vertices,
  and `facets(IncidenceMatrix, P)` for `$p->VERTICES_IN_FACETS`.

## Common polymake commands and their OSCAR counterparts

### Polytopes

| polymake | OSCAR |
|:---------|:------|
| `cube(3)`, `simplex(3)`, `cross(3)` | `cube(3)`, `simplex(3)`, `cross_polytope(3)` |
| `cyclic(3, 6)`, `hypersimplex(2, 4)` | `cyclic_polytope(3, 6)`, `hypersimplex(2, 4)` |
| `rand_sphere(3, 10)` | `rand_spherical_polytope(3, 10)` |
| `johnson_solid(1)`, `dodecahedron()` | `johnson_solid(1)`, `dodecahedron()` |
| `new Polytope(POINTS=>[[1,0,0],[1,1,0]])` | `convex_hull([0 0; 1 0])` |
| `new Polytope(INEQUALITIES=>...)` | `polyhedron(A, b)` |
| `$p->VERTICES`, `$p->FACETS` | `vertices(P)`, `facets(P)` |
| `$p->N_VERTICES`, `$p->N_FACETS` | `n_vertices(P)`, `n_facets(P)` |
| `$p->DIM`, `$p->AMBIENT_DIM` | `dim(P)`, `ambient_dim(P)` |
| `$p->F_VECTOR`, `$p->VOLUME` | `f_vector(P)`, `volume(P)` |
| `$p->SIMPLE`, `$p->SIMPLICIAL` | `is_simple(P)`, `is_simplicial(P)` |
| `$p->BOUNDED`, `$p->FEASIBLE` | `is_bounded(P)`, `is_feasible(P)` |
| `$p->LATTICE`, `$p->LATTICE_POINTS_GENERATORS` | `is_lattice_polytope(P)`, `lattice_points(P)` |
| `$p->INTERIOR_LATTICE_POINTS` | `interior_lattice_points(P)` |
| `$p->EHRHART_POLYNOMIAL`, `$p->H_STAR_VECTOR` | `ehrhart_polynomial(P)`, `h_star_polynomial(P)` |
| `$p->VERTICES_IN_FACETS` | `facets(IncidenceMatrix, P)` |
| `$p->AFFINE_HULL`, `$p->LINEALITY_SPACE` | `affine_hull(P)`, `lineality_space(P)` |
| `$p->GRAPH` | `vertex_edge_graph(P)` |
| `$p->HASSE_DIAGRAM` | `faces(P, k)` |
| `combinatorial_symmetries($p)` | `combinatorial_symmetries(P)` |
| `$p->TRIANGULATION` | `regular_triangulations(P)`, `all_triangulations(P)` |
| `polarize($p)` | `polarize(P)` |
| `minkowski_sum($p, $q)` | `minkowski_sum(P, Q)` or `P + Q` |
| `product($p, $q)` | `product(P, Q)` or `P * Q` |
| `intersection($p, $q)` | `intersect(P, Q)` |
| `scale($p, 2)`, `translate($p, $v)` | `2 * P`, `P + v` |

### Cones and fans

| polymake | OSCAR |
|:---------|:------|
| `new Cone(INPUT_RAYS=>[[1,0],[1,1]])` | `positive_hull([1 0; 1 1])` |
| `$c->RAYS`, `$c->FACETS`, `$c->DIM` | `rays(C)`, `facets(C)`, `dim(C)` |
| `$c->HILBERT_BASIS_GENERATORS` | `hilbert_basis(C)` |
| `normal_fan($p)` (application `fan`) | `normal_fan(P)` |
| `new PolyhedralFan(INPUT_RAYS=>..., INPUT_CONES=>...)` | `polyhedral_fan(IncidenceMatrix(cones), rays)` |
| `$f->RAYS`, `$f->MAXIMAL_CONES` | `rays(F)`, `maximal_cones(F)` |
| `$f->N_RAYS`, `$f->N_MAXIMAL_CONES` | `n_rays(F)`, `n_maximal_cones(F)` |
| `$f->COMPLETE` | `is_complete(F)` |

### Combinatorics

| polymake | OSCAR |
|:---------|:------|
| `uniform_matroid(2, 4)`, `fano_matroid()` (application `matroid`) | `uniform_matroid(2, 4)`, `fano_matroid()` |
| `new Matroid(VECTORS=>...)` | `matroid_from_matrix_columns(M)` |
| `new SimplicialComplex(FACETS=>[[0,1],[1,2]])` (application `topaz`) | `simplicial_complex([[1, 2], [2, 3]])` |
| `graph_from_edges([[0,1],[1,2]])` (application `graph`) | `graph_from_edges([[1, 2], [2, 3]])` |
