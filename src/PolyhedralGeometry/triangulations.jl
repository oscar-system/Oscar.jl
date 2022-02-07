@doc Markdown.doc"""
    all_triangulations(pts::AnyVecOrMat)

Compute all triangulations on the points given as the rows of `pts`.

The return type is a `Vector{Vector{Vector{Int}}}` where each
`Vector{Vector{Int}}` encodes a triangulation, in which a `Vector{Int}` encodes
a simplex as the set of indices of the vertices of the simplex. I.e. the
`Vector{Int}` `[1,2,4]` corresponds to the simplex that is the convex hull of
the first, second, and fourth input point.

# Examples
```jldoctest
julia> c = cube(2,0,1)
A polyhedron in ambient dimension 2

julia> V = vertices(c)
4-element SubObjectIterator{PointVector{Polymake.Rational}}:
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
function all_triangulations(pts::Union{SubObjectIterator{<:PointVector}, AbstractMatrix, Oscar.MatElem})
    input = homogenized_matrix(pts, 1)
    PC = Polymake.polytope.PointConfiguration(POINTS=input)
    triangs = Polymake.polytope.topcom_all_triangulations(PC)
    result = [[[e for e in simplex] for simplex in triang] for triang in triangs]
    return Polymake.to_one_based_indexing(result)
end


@doc Markdown.doc"""
    all_triangulations(P::Polyhedron)

Compute all triangulations that can be formed using the vertices of the given
polytope `P`.

The return type is a `Vector{Vector{Vector{Int}}}` where each
`Vector{Vector{Int}}` encodes a triangulation, in which a `Vector{Int}` encodes
a simplex as the set of indices of the vertices of the simplex. I.e. the
`Vector{Int}` `[1,2,4]` corresponds to the simplex that is the convex hull of
the first, second, and fourth input point.

# Examples
```jldoctest
julia> c = cube(2,0,1)
A polyhedron in ambient dimension 2

julia> all_triangulations(c)
2-element Vector{Vector{Vector{Int64}}}:
 [[1, 2, 3], [2, 3, 4]]
 [[1, 2, 4], [1, 3, 4]]
```
"""
all_triangulations(P::Polyhedron) = all_triangulations(vertices(P))


@doc Markdown.doc"""
    regular_triangulations(pts::AnyVecOrMat)

Compute all regular triangulations on the points given as the rows of `pts`.

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
A polyhedron in ambient dimension 2

julia> V = vertices(c)
4-element SubObjectIterator{PointVector{Polymake.Rational}}:
 [0, 0]
 [1, 0]
 [0, 1]
 [1, 1]

julia> regular_triangulations(V)
2-element Vector{Vector{Vector{Int64}}}:
 [[1, 2, 3], [2, 3, 4]]
 [[1, 3, 4], [1, 2, 4]]
```
"""
function regular_triangulations(pts::Union{SubObjectIterator{<:PointVector}, AbstractMatrix, Oscar.MatElem})
    input = homogenized_matrix(pts, 1)
    PC = Polymake.polytope.PointConfiguration(POINTS=input)
    triangs = Polymake.polytope.topcom_regular_triangulations(PC)
    result = [[[e for e in simplex] for simplex in triang] for triang in triangs]
    return Polymake.to_one_based_indexing(result)
end


@doc Markdown.doc"""
    regular_triangulations(P::Polyhedron)

Compute all regular triangulations that can be formed using the vertices of the
given polytope `P`.

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
A polyhedron in ambient dimension 2

julia> regular_triangulations(c)
2-element Vector{Vector{Vector{Int64}}}:
 [[1, 2, 3], [2, 3, 4]]
 [[1, 3, 4], [1, 2, 4]]
```
"""
regular_triangulations(P::Polyhedron) = regular_triangulations(vertices(P))


@doc Markdown.doc"""
    secondary_polytope(P::Polyhedron)

Compute the secondary polytope of a polyhedron, i.e. the convex hull of all the
gkz vectors of all its (regular) triangulations. A triangulation here means
only using the vertices of `P`.

# Examples
Compute the secondary polytope of the cube.
```jldoctest
julia> c = cube(3)
A polyhedron in ambient dimension 3

julia> sc = secondary_polytope(c)
A polyhedron in ambient dimension 8
```
"""
function secondary_polytope(P::Polyhedron)
    return Polyhedron(Polymake.polytope.secondary_polytope(pm_object(P)))
end
