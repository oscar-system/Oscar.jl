@doc Markdown.doc"""
    ambient_dim(PC::PolyhedralComplex)

Return the ambient dimension of `PC`.

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2,3],[1,3,4]])
2×4 IncidenceMatrix
[1, 2, 3]
[1, 3, 4]


julia> V = [0 0; 1 0; 1 1; 0 1]
4×2 Matrix{Int64}:
 0  0
 1  0
 1  1
 0  1

julia> PC = PolyhedralComplex(IM, V)
A polyhedral complex in ambient dimension 2

julia> ambient_dim(PC)
2
```
"""
ambient_dim(PC::PolyhedralComplex) = Polymake.fan.ambient_dim(pm_object(PC))


@doc Markdown.doc"""
    vertices(as, PC)

Return an iterator over the vertices of `PC` in the format defined by `as`.

Optional arguments for `as` include
* `PointVector`.

# Examples
The following code computes the vertices of the Minkowski sum of a triangle and
a square:
```jldoctest
julia> P = simplex(2) + cube(2);

julia> vertices(PointVector, P)
5-element SubObjectIterator{PointVector{Polymake.Rational}}:
 [-1, -1]
 [2, -1]
 [2, 1]
 [-1, 2]
 [1, 2]
```
"""
vertices(as::Type{PointVector{T}}, PC::PolyhedralComplex) where T = SubObjectIterator{as}(pm_object(PC), _vertex_polyhedron, length(_vertex_indices(pm_object(PC))))


function _all_vertex_indices(P::Polymake.BigObject)
    vi = Polymake.get_attachment(P, "_vertex_indices")
    if isnothing(vi)
        A = P.VERTICES
        vi = Polymake.Vector{Polymake.to_cxx_type(Int64)}(Vector(1:Polymake.nrows(A)))
        Polymake.attach(P, "_all_vertex_indices", vi)
    end
    return vi
end

function _vertex_or_ray_polyhedron(::Type{Union{PointVector{T}, RayVector{T}}}, P::Polymake.BigObject, i::Base.Integer) where T 
    A = P.VERTICES
    if iszero(A[_all_vertex_indices(P)[i],1])
        return RayVector{T}(P.VERTICES[_all_vertex_indices(P)[i], 2:end])
    else
        return PointVector{T}(P.VERTICES[_all_vertex_indices(P)[i], 2:end])
    end
end

vertices(as::Type{Union{RayVector{T}, PointVector{T}}}, PC::PolyhedralComplex) where T = SubObjectIterator{as}(pm_object(PC), _vertex_or_ray_polyhedron, length(_all_vertex_indices(pm_object(PC))))

vertices(::Type{Union{RayVector, PointVector}}, PC::PolyhedralComplex) = vertices(Union{RayVector{Polymake.Rational}, PointVector{Polymake.Rational}}, PC)

vertices(PC::PolyhedralComplex) = vertices(Union{PointVector, RayVector}, PC)

_maximal_polytope(::Type{Polyhedron}, PC::Polymake.BigObject, i::Base.Integer) = Polyhedron(Polymake.fan.polytope(PC, i-1))


@doc Markdown.doc"""
    maximal_polytopes(PC::PolyhedralComplex)

Return the maximal polytopes of `PC`
"""
maximal_polytopes(PC::PolyhedralComplex) = SubObjectIterator{Polyhedron}(pm_object(PC), _maximal_polytope, nmaximal_polytopes(PC))


@doc Markdown.doc"""
    nmaximal_polytopes(PC::PolyhedralComplex)

Return the number of maximal polytopes of `PC`
"""
nmaximal_polytopes(PC::PolyhedralComplex) = pm_object(PC).N_MAXIMAL_POLYTOPES


