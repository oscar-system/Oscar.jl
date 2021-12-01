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


function secondary_polytope(P::Polyhedron)
    return Polymake.polytope.secondary_polytope(pm_object(P))
end
