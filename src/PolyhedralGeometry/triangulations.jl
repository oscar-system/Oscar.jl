################################################################################
# TODO: Remove the following two functions after next polymake release 4.7
function _is_full_triangulation(triang::Vector{Vector{Int}}, npoints::Int)
    u = Set{Int}()
    for v in triang
        union!(u, v)
        if length(u) == npoints
            return true
        end
    end
    return false
end

function _postprocess(triangs::Vector{Vector{Vector{Int}}}, full::Bool)
    result = Polymake.to_one_based_indexing(triangs)
    if full
        result = [t for t in result if _is_full_triangulation(t)]
    end
    return result
end
################################################################################

function _is_star_triangulation(triang::Vector{Vector{Int}})
    u = Set{Int}()
    for v in triang
        if !(1 in v)
            return false
        end
    end
    return true
end

@doc Markdown.doc"""
    all_triangulations(pts::AnyVecOrMat; full=false)

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
function all_triangulations(pts::Union{SubObjectIterator{<:PointVector}, AbstractMatrix, Oscar.MatElem}; full::Bool=false)
    input = homogenized_matrix(pts, 1)
    PC = Polymake.polytope.PointConfiguration(POINTS=input)
    triangs = Polymake.polytope.topcom_all_triangulations(PC)
    result = [[[e for e in simplex] for simplex in triang] for triang in triangs]
    return _postprocess(result, full)
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
all_triangulations(P::Polyhedron) = all_triangulations(vertices(P); full=false)


@doc Markdown.doc"""
    star_triangulations(pts::AnyVecOrMat; full::Bool=false, regular::Bool=false)

Return all star triangulations of the given point configuration, i.e. all
simplices are required to contain the first point. Optionally only select the
`regular` or `full` triangulations.

The return type is a `Vector{Vector{Vector{Int}}}` where each
`Vector{Vector{Int}}` encodes a triangulation, in which a `Vector{Int}` encodes
a simplex as the set of indices of the vertices of the simplex. I.e. the
`Vector{Int}` `[1,2,4]` corresponds to the simplex that is the convex hull of
the first, second, and fourth input point.
"""    
function star_triangulations(pts::AnyVecOrMat; full::Bool=false, regular::Bool=false)
    if regular
        result = regular_triangulations(pts; full=full)
    else
        result = all_triangulations(pts; full=full)
    end
    return [t for t in result if _is_star_triangulation(t)]
end

@doc Markdown.doc"""
    star_triangulations(P::Polyhedron; full::Bool=false, regular::Bool=false)

Return all star triangulations of the given polyhedron, i.e. all simplices are
required to contain the origin. If the origin is not among the vertices of the
polyhedron, it is added. Optionally only select the `regular` or `full`
triangulations.

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
A polyhedron in ambient dimension 2

julia> star_triangulations(hex)
([0 0; -1 -1; 0 -1; 1 0; 1 1; 0 1; -1 0], [[[1, 2, 3], [1, 2, 7], [1, 3, 4], [1, 4, 5], [1, 5, 6], [1, 6, 7]]])
```
"""    
function star_triangulations(P::Polyhedron; full::Bool=false, regular::Bool=false)
    zero = [0 for i in 1:dim(P)]
    contains(P, zero) || throw(ArgumentError("Input polyhedron must contain origin."))
    V = vertices(P)
    if zero in V
        V = [v for v in V if v != 0]
    end
    pts = vcat(matrix(QQ, transpose(zero)), matrix(QQ, V))
    return pts, star_triangulations(pts; full=full, regular=regular)
end




@doc Markdown.doc"""
    regular_triangulations(pts::AnyVecOrMat; full=false)

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
function regular_triangulations(pts::Union{SubObjectIterator{<:PointVector}, AbstractMatrix, Oscar.MatElem}; full::Bool=false)
    input = homogenized_matrix(pts, 1)
    PC = Polymake.polytope.PointConfiguration(POINTS=input)
    triangs = Polymake.polytope.topcom_regular_triangulations(PC)
    result = [[[e for e in simplex] for simplex in triang] for triang in triangs]
    return _postprocess(result, full)
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
regular_triangulations(P::Polyhedron) = regular_triangulations(vertices(P); full=false)


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
