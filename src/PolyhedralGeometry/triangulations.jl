import JSON
using TOPCOM_jll

function _postprocess_polymake_triangs(
  triangs::Polymake.Array{Polymake.Set{Polymake.Set{Polymake.to_cxx_type(Int)}}}
)
  result = Vector{Vector{Int}}[]
  for triang in triangs
    push!(result, [Polymake.to_one_based_indexing(Vector{Int}(t)) for t in triang])
  end
  return result
end

function topcom_regular_triangulations(
  pts::AbstractCollection[PointVector]; full::Bool=false
)
  input = homogenized_matrix(pts, 1)
  PC = Polymake.polytope.PointConfiguration(; POINTS=input)
  result = if full
    Polymake.polytope.topcom_fine_and_regular_triangulations(PC)
  else
    Polymake.polytope.topcom_regular_triangulations(PC)
  end
  return _postprocess_polymake_triangs(result)
end

function topcom_regular_triangulation(
  pts::AbstractCollection[PointVector]; full::Bool=false
)
  input = homogenized_matrix(pts, 1)
  inputstr = join(["[" * join(input[i, :], ",") * "]" for i in 1:nrows(input)], ",\n")
  in = Pipe()
  out = Pipe()
  err = Pipe()
  Base.link_pipe!(in; writer_supports_async=true)
  Base.link_pipe!(out; reader_supports_async=true)
  Base.link_pipe!(err; reader_supports_async=true)
  cmd = Oscar.TOPCOM_jll.points2placingtriang()
  if full
    cmd = Oscar.TOPCOM_jll.points2finetriang()
  end
  proc = run(pipeline(`$(cmd) --regular`; stdin=in, stdout=out, stderr=err); wait=false)
  task = @async begin
    write(in, "[\n$inputstr\n]\n")
    close(in)
  end
  close(in.out)
  close(out.in)
  close(err.in)
  result = Vector{Vector{Vector{Int}}}()
  for line in eachline(out)
    m = match(r"{{.*}}", line)
    triang = replace(m.match, "{" => "[")
    triang = replace(triang, "}" => "]")
    triang = convert(Vector{Vector{Int}}, JSON.parse(triang))
    push!(result, Polymake.to_one_based_indexing(triang))
  end
  wait(task)
  if !success(proc)
    msg = eof(err) ? "unknown error" : readchomp(err)
    error("Failed to run TOPCOM: $msg")
  end
  return result
end

################################################################################

function _is_triangulation(sop::SubdivisionOfPoints{QQFieldElem})
  ad = ambient_dim(sop)
  for mc in maximal_cells(sop)
    ad == length(mc) - 1 || return false
  end
  return true
end

function _is_full_triangulation(sop::SubdivisionOfPoints{QQFieldElem})
  _is_triangulation(sop) || return false
  length(Base.union(maximal_cells(sop)...)) == n_points(sop) || return false
  return true
end

_is_star_triangulation(sop::SubdivisionOfPoints{QQFieldElem}) =
  all(cell -> 1 in cell, maximal_cells(sop))

function _find_full_star_triangulation(pts::ZZMatrix; seed::Int=-1)
  seed == -1 || Random.seed!(seed)
  n = nrows(pts)
  wts = rand(0:100000, n)
  # Weight of first points is lowest
  wts[1] = -1000000
  sop = subdivision_of_points(pts, wts)
  while !(_is_full_triangulation(sop) && _is_star_triangulation(sop))
    wts = rand(0:100000, n)
    wts[1] = -1000000
    sop = subdivision_of_points(pts, wts)
  end
  return collect(maximal_cells(sop))
end

function _is_star_triangulation(triang::Vector{Vector{Int}})
  u = Set{Int}()
  for v in triang
    if !(1 in v)
      return false
    end
  end
  return true
end

@doc raw"""
    all_triangulations(pts::AbstractCollection[PointVector]; full=false)

Compute all triangulations on the points given as the rows of `pts`. Optionally
select `full=true` to output full triangulations only, i.e. those that use all
given points.

The return type is a `Vector{Vector{Vector{Int}}}` where each
`Vector{Vector{Int}}` encodes a triangulation, in which a `Vector{Int}` encodes
a simplex as the set of indices of the vertices of the simplex. I.e. the
`Vector{Int}` `[1,2,4]` corresponds to the simplex that is the convex hull of
the first, second, and fourth input point.

# Examples
```jldoctest
julia> c = cube(2,0,1)
Polytope in ambient dimension 2

julia> V = vertices(c)
4-element SubObjectIterator{PointVector{QQFieldElem}}:
 [0, 0]
 [1, 0]
 [0, 1]
 [1, 1]

julia> all_triangulations(V)
2-element Vector{Vector{Vector{Int64}}}:
 [[1, 2, 3], [2, 3, 4]]
 [[1, 2, 4], [1, 3, 4]]
```
"""
function all_triangulations(pts::AbstractCollection[PointVector]; full::Bool=false)
  input = homogenized_matrix(pts, 1)
  PC = Polymake.polytope.PointConfiguration(; POINTS=input)
  PC.FULL_DIM::Bool || error("Input points must have full rank.")
  result = if full
    Polymake.polytope.topcom_fine_triangulations(PC)
  else
    Polymake.polytope.topcom_all_triangulations(PC)
  end
  return _postprocess_polymake_triangs(result)
end

@doc raw"""
    all_triangulations(P::Polyhedron)

Compute all triangulations that can be formed using the vertices of the given
bounded and full-dimensional polytope `P`.

The return type is a `Vector{Vector{Vector{Int}}}` where each
`Vector{Vector{Int}}` encodes a triangulation, in which a `Vector{Int}` encodes
a simplex as the set of indices of the vertices of the simplex. I.e. the
`Vector{Int}` `[1,2,4]` corresponds to the simplex that is the convex hull of
the first, second, and fourth input point.

# Examples
```jldoctest
julia> c = cube(2,0,1)
Polytope in ambient dimension 2

julia> all_triangulations(c)
2-element Vector{Vector{Vector{Int64}}}:
 [[1, 2, 3], [2, 3, 4]]
 [[1, 2, 4], [1, 3, 4]]
```
"""
function all_triangulations(P::Polyhedron)
  is_fulldimensional(P) || error("Input polytope must be full-dimensional.")
  is_bounded(P) || error("Input polytope must be bounded.")
  return all_triangulations(vertices(P); full=false)
end

@doc raw"""
    star_triangulations(pts::AbstractCollection[PointVector]; full::Bool=false, regular::Bool=false)

Return all star triangulations of the given point configuration, i.e. all
simplices are required to contain the first point. Optionally only select the
`regular` or `full` triangulations.

The return type is a `Vector{Vector{Vector{Int}}}` where each
`Vector{Vector{Int}}` encodes a triangulation, in which a `Vector{Int}` encodes
a simplex as the set of indices of the vertices of the simplex. I.e. the
`Vector{Int}` `[1,2,4]` corresponds to the simplex that is the convex hull of
the first, second, and fourth input point.
"""
function star_triangulations(
  pts::AbstractCollection[PointVector]; full::Bool=false, regular::Bool=false
)
  if regular
    result = regular_triangulations(pts; full=full)
  else
    result = all_triangulations(pts; full=full)
  end
  return [t for t in result if _is_star_triangulation(t)]
end

@doc raw"""
    star_triangulations(P::Polyhedron; full::Bool=false, regular::Bool=false)

Return all star triangulations of the given bounded and full-dimensional
polyhedron, i.e. all simplices are required to contain the origin. If the
origin is not among the vertices of the polyhedron, it is added. Optionally
only select the `regular` or `full` triangulations.

The output is a pair
- The first entry contains the points of the point configuration as their order
may be different than inside the polyhedron.
- The second entry contains the triangulations as
`Vector{Vector{Vector{Int}}}`, where each `Vector{Vector{Int}}` encodes a
triangulation, in which a `Vector{Int}` encodes a simplex as the set of indices
of the vertices of the simplex. I.e. the `Vector{Int}` `[1,2,4]` corresponds to
the simplex that is the convex hull of the first, second, and fourth input
point.

# Examples
A two-dimensional polyhedron has only one star triangulation.
```jldoctest
julia> hex = convex_hull([-1 -1; 0 -1; 1 0; 1 1; 0 1; -1 0])
Polyhedron in ambient dimension 2

julia> star_triangulations(hex)
([0 0; -1 -1; 0 -1; 1 0; 1 1; 0 1; -1 0], [[[1, 2, 3], [1, 2, 7], [1, 3, 4], [1, 4, 5], [1, 5, 6], [1, 6, 7]]])

julia> star_triangulations(hex; full=true)
([0 0; -1 -1; 0 -1; 1 0; 1 1; 0 1; -1 0], [[[1, 2, 3], [1, 2, 7], [1, 3, 4], [1, 4, 5], [1, 5, 6], [1, 6, 7]]])

julia> star_triangulations(hex; full=true, regular=true)
([0 0; -1 -1; 0 -1; 1 0; 1 1; 0 1; -1 0], [[[1, 2, 3], [1, 2, 7], [1, 3, 4], [1, 4, 5], [1, 5, 6], [1, 6, 7]]])
```
A three-dimensional example with two star triangulations.
```jldoctest
julia> P = convex_hull([0 0 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1])
Polyhedron in ambient dimension 3

julia> star_triangulations(P)
([0 0 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1], [[[1, 2, 3, 4], [1, 2, 4, 5]], [[1, 2, 3, 5], [1, 3, 4, 5]]])
```
"""
function star_triangulations(P::Polyhedron; full::Bool=false, regular::Bool=false)
  is_fulldimensional(P) || error("Input polytope must be full-dimensional.")
  is_bounded(P) || error("Input polytope must be bounded.")
  zero = [0 for i in 1:ambient_dim(P)]
  @req zero in P "Input polyhedron must contain origin."
  V = vertices(P)
  V = [Vector{QQFieldElem}(v) for v in V if !iszero(v)]
  pts = vcat(matrix(QQ, transpose(zero)), matrix(QQ, transpose(hcat(V...))))
  return pts, star_triangulations(pts; full=full, regular=regular)
end

@doc raw"""
    regular_triangulations(pts::AbstractCollection[PointVector]; full=false)

Compute all regular triangulations on the points given as the rows of `pts`.

A triangulation is regular if it can be induced by weights, i.e. attach a
weight to every point, take the convex hull of these new vectors and then take
the subdivision corresponding to the facets visible from below (lower
envelope). Optionally specify `full`, i.e. that every triangulation must use
all points.

The return type is a `Vector{Vector{Vector{Int}}}` where each
`Vector{Vector{Int}}` encodes a triangulation, in which a `Vector{Int}` encodes
a simplex as the set of indices of the vertices of the simplex. I.e. the
`Vector{Int}` `[1,2,4]` corresponds to the simplex that is the convex hull of
the first, second, and fourth input point.

# Examples
```jldoctest
julia> c = cube(2,0,1)
Polytope in ambient dimension 2

julia> V = vertices(c)
4-element SubObjectIterator{PointVector{QQFieldElem}}:
 [0, 0]
 [1, 0]
 [0, 1]
 [1, 1]

julia> regular_triangulations(V)
2-element Vector{Vector{Vector{Int64}}}:
 [[1, 2, 3], [2, 3, 4]]
 [[1, 2, 4], [1, 3, 4]]
```
"""
function regular_triangulations(pts::AbstractCollection[PointVector]; full::Bool=false)
  input = homogenized_matrix(pts, 1)
  PC = Polymake.polytope.PointConfiguration(; POINTS=input)
  PC.FULL_DIM::Bool || error("Input points must have full rank.")
  return topcom_regular_triangulations(pts; full=full)
end

@doc raw"""
    regular_triangulations(P::Polyhedron)

Compute all regular triangulations that can be formed using the vertices of the
given bounded and full-dimensional polytope `P`.

A triangulation is regular if it can be induced by weights, i.e. attach a
weight to every point, take the convex hull of these new vectors and then take
the subdivision corresponding to the facets visible from below (lower
envelope).

The return type is a `Vector{Vector{Vector{Int}}}` where each
`Vector{Vector{Int}}` encodes a triangulation, in which a `Vector{Int}` encodes
a simplex as the set of indices of the vertices of the simplex. I.e. the
`Vector{Int}` `[1,2,4]` corresponds to the simplex that is the convex hull of
the first, second, and fourth input point.

# Examples
```jldoctest
julia> c = cube(2,0,1)
Polytope in ambient dimension 2

julia> regular_triangulations(c)
2-element Vector{Vector{Vector{Int64}}}:
 [[1, 2, 3], [2, 3, 4]]
 [[1, 2, 4], [1, 3, 4]]
```
"""
function regular_triangulations(P::Polyhedron)
  is_fulldimensional(P) || error("Input polytope must be full-dimensional.")
  is_bounded(P) || error("Input polytope must be bounded.")
  return regular_triangulations(vertices(P); full=false)
end

@doc raw"""
    regular_triangulation(pts::AbstractCollection[PointVector]; full=false)

Compute ONE regular triangulation on the points given as the rows of `pts`.

A triangulation is regular if it can be induced by weights, i.e. attach a
weight to every point, take the convex hull of these new vectors and then take
the subdivision corresponding to the facets visible from below (lower
envelope). Optionally specify `full`, i.e. that every triangulation must use
all points.

As for `regular_triangulation(pts::AnyVecOrMat; full=false)` the return type is
`Vector{Vector{Vector{Int}}}`. Here, only one triangulation is computed, so
the outer vector is of length one. Its entry of type `Vector{Vector{Int}}`
encodes the triangulation in question. Recall that a `Vector{Int}` encodes
a simplex as the set of indices of the vertices of the simplex. I.e. the
`Vector{Int}` `[1,2,4]` corresponds to the simplex that is the convex hull of
the first, second, and fourth input point.

# Examples
```jldoctest
julia> c = cube(2,0,1)
Polytope in ambient dimension 2

julia> V = vertices(c)
4-element SubObjectIterator{PointVector{QQFieldElem}}:
 [0, 0]
 [1, 0]
 [0, 1]
 [1, 1]

julia> regular_triangulation(V)
1-element Vector{Vector{Vector{Int64}}}:
 [[1, 2, 3], [2, 3, 4]]
```
"""
function regular_triangulation(pts::AbstractCollection[PointVector]; full::Bool=false)
  input = homogenized_matrix(pts, 1)
  PC = Polymake.polytope.PointConfiguration(; POINTS=input)
  PC.FULL_DIM::Bool || error("Input points must have full rank.")
  return topcom_regular_triangulation(pts; full=full)
end

@doc raw"""
    regular_triangulation(P::Polyhedron)

Compute ONE regular triangulation that can be formed using the vertices of the
given bounded and full-dimensional polytope `P`.

A triangulation is regular if it can be induced by weights, i.e. attach a
weight to every point, take the convex hull of these new vectors and then take
the subdivision corresponding to the facets visible from below (lower
envelope).

As for `regular_triangulations(P::Polyhedron)` the return type is
`Vector{Vector{Vector{Int}}}`. Here, only one triangulation is computed, so
the outer vector is of length one. Its entry of type `Vector{Vector{Int}}`
encodes the triangulation in question. Recall that a `Vector{Int}` encodes
a simplex as the set of indices of the vertices of the simplex. I.e. the
`Vector{Int}` `[1,2,4]` corresponds to the simplex that is the convex hull of
the first, second, and fourth input point.

# Examples
```jldoctest
julia> c = cube(2,0,1)
Polytope in ambient dimension 2

julia> regular_triangulation(c)
1-element Vector{Vector{Vector{Int64}}}:
 [[1, 2, 3], [2, 3, 4]]
```
"""
function regular_triangulation(P::Polyhedron)
  is_fulldimensional(P) || error("Input polytope must be full-dimensional.")
  is_bounded(P) || error("Input polytope must be bounded.")
  return regular_triangulation(vertices(P); full=false)
end

@doc raw"""
    secondary_polytope(P::Polyhedron)

Compute the secondary polytope of a polyhedron, i.e. the convex hull of all the
gkz vectors of all its (regular) triangulations. A triangulation here means
only using the vertices of `P`.

# Examples
Compute the secondary polytope of the cube.
```jldoctest
julia> c = cube(3)
Polytope in ambient dimension 3

julia> sc = secondary_polytope(c)
Polytope in ambient dimension 8
```
"""
function secondary_polytope(P::Polyhedron{T}) where {T<:scalar_types}
  @req is_bounded(P) "secondary_polytope only works for polytopes, i.e. bounded polyhedra."
  @req is_fulldimensional(P) "Polytope must be full-dimensional, use a projection."
  return Polyhedron{T}(
    Polymake.polytope.secondary_polytope(pm_object(P)), coefficient_field(P)
  )
end

@doc raw"""
    is_regular(pts::AbstractCollection[PointVector], cells::Vector{Vector{Vector{Int64}}})

Compute whether a triangulation is regular.

# Examples
Compute whether a triangulation of the square is regular.
```jldoctest
julia> c = cube(2)
Polytope in ambient dimension 2

julia> cells=[[1,2,3],[2,3,4]];

julia> is_regular(vertices(c),cells)
true
```
"""
function is_regular(pts::AbstractCollection[PointVector], cells::Vector{Vector{Int64}})
  as_sop = subdivision_of_points(pts, cells)
  is_regular(as_sop)
end

@doc raw"""
    subdivision_of_points(P::Polyhdron, cells::IncidenceMatrix)

# Arguments
- `P::Polyhedron`: A polyhedron whose vertices are the points of the subdivision.
- `cells::IncidenceMatrix`: An incidence matrix; there is a 1 at position (i,j) if cell i contains point j, and 0 otherwise.

A subdivision of points formed from points and cells made of these points. The
cells are given as an IncidenceMatrix, where the columns represent the points
and the rows represent the cells.

# Examples
Compute a triangulation of the square
```jldoctest
julia> C = cube(2);

julia> cells = incidence_matrix([[1,2,3],[2,3,4]]);

julia> S = subdivision_of_points(C, cells)
Subdivision of points in ambient dimension 2
```
"""
subdivision_of_points(P::Polyhedron, cells::IncidenceMatrix) =
  subdivision_of_points(vertices(P), cells)

@doc raw"""
    subdivision_of_points(P::Polyhdron, weights::AbstractVector)

# Arguments
- `P::Polyhedron`: A polyhedron whose vertices are the points of the subdivision.
- `weights::AbstractVector`: A vector with one entry for every point indicating the height of this point.

A subdivision of points formed by placing every vertex of `P` at the corresponding
height, then taking the convex hull and then only considering those cells
corresponding to faces visible from below ("lower envelope").

# Examples
Compute a triangulation of the square
```jldoctest
julia> C = cube(2);

julia> weights = [0,0,1,2];

julia> S = subdivision_of_points(C, weights)
Subdivision of points in ambient dimension 2
```
"""
subdivision_of_points(P::Polyhedron, weights::AbstractVector) =
  subdivision_of_points(vertices(P), weights)
subdivision_of_points(P::Polyhedron, cells::Vector{Vector{Int64}}) =
  subdivision_of_points(vertices(P), IncidenceMatrix(cells))
subdivision_of_points(Iter::SubObjectIterator{<:PointVector}, cells::IncidenceMatrix) =
  subdivision_of_points(point_matrix(Iter), cells)
subdivision_of_points(Iter::SubObjectIterator{<:PointVector}, weights::AbstractVector) =
  subdivision_of_points(point_matrix(Iter), weights)
subdivision_of_points(
  Iter::SubObjectIterator{<:PointVector}, cells::Vector{Vector{Int64}}
) = subdivision_of_points(point_matrix(Iter), IncidenceMatrix(cells))

@doc raw"""
    gkz_vector(SOP::SubdivisionOfPoints)

Compute the gkz vector of a triangulation given as a subdivision of points, SOP.

# Examples
Compute the gkz vector of one of the two regular triangulations of the square.
```jldoctest
julia> C = cube(2);

julia> Triang = subdivision_of_points(C,[[1,2,3],[2,3,4]])
Subdivision of points in ambient dimension 2

julia> gkz_vector(Triang)
4-element Vector{QQFieldElem}:
 4
 8
 8
 4
```
"""
function gkz_vector(SOP::SubdivisionOfPoints)
  V = SOP.pm_subdivision.POINTS
  T = SOP.pm_subdivision.MAXIMAL_CELLS
  n = ambient_dim(SOP)
  for i in 1:size(T, 1)
    @assert sum(T[i, :]) == n + 1 #poor check that subdivision is triangulation
  end
  TT = [Polymake.to_zero_based_indexing(Polymake.row(T, i)) for i in 1:Polymake.nrows(T)]
  F = coefficient_field(SOP)
  tt = elem_type(F)
  tt[F(e) for e in Polymake.call_function(:polytope, :gkz_vector, V, TT)]
end
