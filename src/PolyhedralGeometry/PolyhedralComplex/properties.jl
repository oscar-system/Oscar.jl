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

_maximal_polyhedron(::Type{Polyhedron}, PC::Polymake.BigObject, i::Base.Integer) = Polyhedron(Polymake.fan.polytope(PC, i-1))


@doc Markdown.doc"""
    maximal_polyhedra(PC::PolyhedralComplex)

Return the maximal polyhedra of `PC`

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2,3],[1,3,4]])
2×4 IncidenceMatrix
[1, 2, 3]
[1, 3, 4]


julia> VR = [0 0; 1 0; 1 1; 0 1]
4×2 Matrix{Int64}:
 0  0
 1  0
 1  1
 0  1

julia> PC = PolyhedralComplex(IM, VR, [2])
A polyhedral complex in ambient dimension 2

julia> maximal_polyhedra(PC)
2-element SubObjectIterator{Polyhedron}:
 A polyhedron in ambient dimension 2
 A polyhedron in ambient dimension 2
```
"""
maximal_polyhedra(PC::PolyhedralComplex) = SubObjectIterator{Polyhedron}(pm_object(PC), _maximal_polyhedron, nmaximal_polyhedra(PC))


@doc Markdown.doc"""
    nmaximal_polyhedra(PC::PolyhedralComplex)

Return the number of maximal polyhedra of `PC`

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2,3],[1,3,4]])
2×4 IncidenceMatrix
[1, 2, 3]
[1, 3, 4]


julia> VR = [0 0; 1 0; 1 1; 0 1]
4×2 Matrix{Int64}:
 0  0
 1  0
 1  1
 0  1

julia> PC = PolyhedralComplex(IM, VR, [2])
A polyhedral complex in ambient dimension 2

julia> nmaximal_polyhedra(PC)
2
```
"""
nmaximal_polyhedra(PC::PolyhedralComplex) = pm_object(PC).N_MAXIMAL_POLYTOPES


@doc Markdown.doc"""
    issimplicial(PC::PolyhedralComplex)

Determine whether the polyhedral complex is simplicial.

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2,3],[1,3,4]]);

julia> VR = [0 0; 1 0; 1 1; 0 1];

julia> PC = PolyhedralComplex(IM, VR)
A polyhedral complex in ambient dimension 2

julia> issimplicial(PC)
true
```
"""
issimplicial(PC::PolyhedralComplex) = pm_object(PC).SIMPLICIAL::Bool


@doc Markdown.doc"""
    ispure(PC::PolyhedralComplex)

Determine whether the polyhedral complex is pure.

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2,3],[1,3,4]]);

julia> VR = [0 0; 1 0; 1 1; 0 1];

julia> PC = PolyhedralComplex(IM, VR)
A polyhedral complex in ambient dimension 2

julia> ispure(PC)
true
```
"""
ispure(PC::PolyhedralComplex) = pm_object(PC).PURE::Bool


@doc Markdown.doc"""
    dim(PC::PolyhedralComplex)

Compute the dimension of the polyhedral complex.

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2,3],[1,3,4]]);

julia> VR = [0 0; 1 0; 1 1; 0 1];

julia> PC = PolyhedralComplex(IM, VR)
A polyhedral complex in ambient dimension 2

julia> dim(PC)
2
```
"""
dim(PC::PolyhedralComplex) = Polymake.fan.dim(pm_object(PC))

@doc Markdown.doc"""
    polyhedra_of_dim(PC::PolyhedralComplex, polyhedron_dim::Int)

Return the polyhedra of a given dimension in the polyhedral complex `PC`.

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2,3],[1,3,4]]);

julia> VR = [0 0; 1 0; 1 1; 0 1];

julia> PC = PolyhedralComplex(IM, VR);

julia> P1s = polyhedra_of_dim(PC,1)
5-element SubObjectIterator{Polyhedron}:
 A polyhedron in ambient dimension 2
 A polyhedron in ambient dimension 2
 A polyhedron in ambient dimension 2
 A polyhedron in ambient dimension 2
 A polyhedron in ambient dimension 2

julia> for p in P1s
       println(dim(p))
       end
1
1
1
1
1
```
"""
function polyhedra_of_dim(PC::PolyhedralComplex, polyhedron_dim::Int)
    n = polyhedron_dim - lineality_dim(PC) + 1
    n < 0 && return nothing
    pfaces = Polymake.fan.cones_of_dim(pm_object(PC), n)
    nfaces = Polymake.nrows(pfaces)
    rfaces = Vector{Int64}()
    nfarf = 0
    farf = Polymake.to_one_based_indexing(pm_object(PC).FAR_VERTICES)
    for index in 1:nfaces
        face = Polymake.row(pfaces, index)
        if face <= farf
            nfarf += 1
        else
            append!(rfaces, index)
        end
    end
    return SubObjectIterator{Polyhedron}(pm_object(PC), _ith_polyhedron, length(rfaces), (f_dim = n, f_ind = rfaces))
end

function _ith_polyhedron(::Type{Polyhedron}, PC::Polymake.BigObject, i::Base.Integer; f_dim::Int = -1, f_ind::Vector{Int64} = Vector{Int64}())
    pface = Polymake.row(Polymake.fan.cones_of_dim(PC, f_dim), f_ind[i])
    return Polyhedron(Polymake.polytope.Polytope(VERTICES = PC.VERTICES[collect(pface),:], LINEALITY_SPACE = PC.LINEALITY_SPACE))
end

lineality_space(PC::PolyhedralComplex) = SubObjectIterator{RayVector{Polymake.Rational}}(pm_object(PC), _lineality_polyhedron, lineality_dim(PC))

lineality_dim(PC::PolyhedralComplex) = pm_object(PC).LINEALITY_DIM
