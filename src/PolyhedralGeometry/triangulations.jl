@doc Markdown.doc"""
    all_triangulations(pts::AnyVecOrMat)

Compute all triangulations on the points given as the rows of `pts`.

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

julia> all_triangulations(point_matrix(V))
pm::Array<pm::Array<pm::Set<long, pm::operations::cmp>>>
<{0 1 2}
{1 2 3}
>
<{0 1 3}
{0 2 3}
>
```
"""
function all_triangulations(pts::AnyVecOrMat)
    input = homogenize(matrix_for_polymake(pts))
    PC = Polymake.polytope.PointConfiguration(POINTS=input)
    return Polymake.polytope.topcom_all_triangulations(PC)
end


@doc Markdown.doc"""
    all_triangulations(P::Polyhedron)

Compute all triangulations that can be formed using the vertices of the given
polytope `P`.

# Examples
```jldoctest
julia> c = cube(2,0,1)
A polyhedron in ambient dimension 2

julia> all_triangulations(c)
pm::Array<pm::Array<pm::Set<long, pm::operations::cmp>>>
<{0 1 2}
{1 2 3}
>
<{0 1 3}
{0 2 3}
>
```
"""
all_triangulations(P::Polyhedron) = all_triangulations(point_matrix(vertices(P)))


@doc Markdown.doc"""
    regular_triangulations(pts::AnyVecOrMat)

Compute all regular triangulations on the points given as the rows of `pts`.

A triangulation is regular if it can be induced by weights, i.e. attach a
weight to every point, take the convex hull of these new vectors and then take
the subdivision corresponding to the facets visible from below (lower
envelope).

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

julia> regular_triangulations(point_matrix(V))
pm::Array<pm::Array<pm::Set<long, pm::operations::cmp>>>
<{0 1 2}
{1 2 3}
>
<{0 1 3}
{0 2 3}
>
```
"""
function regular_triangulations(pts::AnyVecOrMat)
    input = homogenize(matrix_for_polymake(pts))
    PC = Polymake.polytope.PointConfiguration(POINTS=input)
    return Polymake.polytope.topcom_regular_triangulations(PC)
end


@doc Markdown.doc"""
    regular_triangulations(P::Polyhedron)

Compute all regular triangulations that can be formed using the vertices of the
given polytope `P`.

A triangulation is regular if it can be induced by weights, i.e. attach a
weight to every point, take the convex hull of these new vectors and then take
the subdivision corresponding to the facets visible from below (lower
envelope).

# Examples
```jldoctest
julia> c = cube(2,0,1)
A polyhedron in ambient dimension 2

julia> regular_triangulations(c)
pm::Array<pm::Array<pm::Set<long, pm::operations::cmp>>>
<{0 1 2}
{1 2 3}
>
<{0 1 3}
{0 2 3}
>
```
"""
regular_triangulations(P::Polyhedron) = regular_triangulations(point_matrix(vertices(P)))


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
