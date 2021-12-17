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
ambient_dim(PC::PolyhedralComplex) = Polymake.fan.ambient_dim(pm_object(PC))::Int


@doc Markdown.doc"""
    vertices([as::Type,] PC::PolyhedralComplex)

Return an iterator over the vertices of `PC` in the format defined by `as`.

Optional arguments for `as` include
* `PointVector`.

"""
vertices(as::Type{PointVector{T}}, PC::PolyhedralComplex) where T<:scalar_types = SubObjectIterator{as}(pm_object(PC), _vertex_polyhedron, length(_vertex_indices(pm_object(PC))))


function _all_vertex_indices(P::Polymake.BigObject)
    vi = Polymake.get_attachment(P, "_all_vertex_indices")
    if isnothing(vi)
        A = P.VERTICES
        vi = Polymake.Vector{Polymake.to_cxx_type(Int64)}(Vector(1:Polymake.nrows(A)))
        Polymake.attach(P, "_all_vertex_indices", vi)
    end
    return vi
end

function _vertex_or_ray_polyhedron(::Type{Union{PointVector{T}, RayVector{T}}}, P::Polymake.BigObject, i::Base.Integer) where T<:scalar_types
    A = P.VERTICES
    if iszero(A[_all_vertex_indices(P)[i],1])
        return RayVector{T}(P.VERTICES[_all_vertex_indices(P)[i], 2:end])
    else
        return PointVector{T}(P.VERTICES[_all_vertex_indices(P)[i], 2:end])
    end
end

vertices(as::Type{Union{RayVector{T}, PointVector{T}}}, PC::PolyhedralComplex) where T<:scalar_types = SubObjectIterator{as}(pm_object(PC), _vertex_or_ray_polyhedron, length(_all_vertex_indices(pm_object(PC))))

vertices(::Type{Union{RayVector, PointVector}}, PC::PolyhedralComplex) where T<:scalar_types = vertices(Union{RayVector{T}, PointVector{T}}, PC)

vertices_and_rays(PC::PolyhedralComplex{T}) where T<:scalar_types = vertices(Union{PointVector{T}, RayVector{T}}, PC)

_vector_matrix(::Val{_vertex_or_ray_polyhedron}, PC::Polymake.BigObject; homogenized = false) = PC.VERTICES[:, (homogenized ? 1 : 2):end]

vertices(::Type{PointVector}, PC::PolyhedralComplex{T}) where T<:scalar_types = vertices(PointVector{T}, PC)

vertices(PC::PolyhedralComplex) = vertices(PointVector, PC)

rays(as::Type{RayVector{T}}, PC::PolyhedralComplex) where T<:scalar_types = SubObjectIterator{as}(pm_object(PC), _ray_polyhedral_complex, length(_ray_indices_polyhedral_complex(pm_object(PC))))

rays(::Type{RayVector}, PC::PolyhedralComplex{T}) where T<:scalar_types = rays(RayVector{T}, PC)

@doc Markdown.doc"""
    rays(PC::PolyhedralComplex)

Return the rays of `PC`

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2,3],[1,3,4]]);

julia> VR = [0 0; 1 0; 1 1; 0 1];

julia> PC = PolyhedralComplex(IM, VR, [2])
A polyhedral complex in ambient dimension 2

julia> rays(PC);

julia> rays(PC)
1-element SubObjectIterator{RayVector{fmpq}}:
 [1, 0]
```
"""
rays(PC::PolyhedralComplex) = rays(RayVector,PC)

_ray_indices_polyhedral_complex(PC::Polymake.BigObject) = collect(Polymake.to_one_based_indexing(PC.FAR_VERTICES))

_ray_polyhedral_complex(::Type{RayVector{T}}, PC::Polymake.BigObject, i::Base.Integer) where T<:scalar_types = RayVector{T}(PC.VERTICES[_ray_indices_polyhedral_complex(PC)[i], 2:end])

_vector_matrix(::Val{_ray_polyhedral_complex}, PC::Polymake.BigObject; homogenized = false) where T<:scalar_types = PC.VERTICES[_ray_indices_polyhedral_complex(PC), (homogenized ? 1 : 2):end]

_maximal_polyhedron(::Type{Polyhedron{T}}, PC::Polymake.BigObject, i::Base.Integer) where T<:scalar_types = Polyhedron{T}(Polymake.fan.polytope(PC, i-1))


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
2-element SubObjectIterator{Polyhedron{fmpq}}:
 A polyhedron in ambient dimension 2
 A polyhedron in ambient dimension 2
```
"""
maximal_polyhedra(PC::PolyhedralComplex{T}) where T<:scalar_types = SubObjectIterator{Polyhedron{T}}(pm_object(PC), _maximal_polyhedron, n_maximal_polyhedra(PC))


@doc Markdown.doc"""
    n_maximal_polyhedra(PC::PolyhedralComplex)

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

julia> n_maximal_polyhedra(PC)
2
```
"""
n_maximal_polyhedra(PC::PolyhedralComplex) = pm_object(PC).N_MAXIMAL_POLYTOPES


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
dim(PC::PolyhedralComplex) = Polymake.fan.dim(pm_object(PC))::Int

@doc Markdown.doc"""
    polyhedra_of_dim(PC::PolyhedralComplex, polyhedron_dim::Int)

Return the polyhedra of a given dimension in the polyhedral complex `PC`.

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2,3],[1,3,4]]);

julia> VR = [0 0; 1 0; 1 1; 0 1];

julia> PC = PolyhedralComplex(IM, VR);

julia> P1s = polyhedra_of_dim(PC,1)
5-element SubObjectIterator{Polyhedron{fmpq}}:
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
function polyhedra_of_dim(PC::PolyhedralComplex{T}, polyhedron_dim::Int) where T<:scalar_types
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
    return SubObjectIterator{Polyhedron{T}}(pm_object(PC), _ith_polyhedron, length(rfaces), (f_dim = n, f_ind = rfaces))
end

function _ith_polyhedron(::Type{Polyhedron{T}}, PC::Polymake.BigObject, i::Base.Integer; f_dim::Int = -1, f_ind::Vector{Int64} = Vector{Int64}()) where T<:scalar_types
    pface = Polymake.row(Polymake.fan.cones_of_dim(PC, f_dim), f_ind[i])
    return Polyhedron{T}(Polymake.polytope.Polytope{scalar_type_to_polymake[T]}(VERTICES = PC.VERTICES[collect(pface),:], LINEALITY_SPACE = PC.LINEALITY_SPACE))
end

lineality_space(PC::PolyhedralComplex{T}) where T<:scalar_types = SubObjectIterator{RayVector{T}}(pm_object(PC), _lineality_polyhedron, lineality_dim(PC))

lineality_dim(PC::PolyhedralComplex) = pm_object(PC).LINEALITY_DIM::Int


@doc Markdown.doc"""
    f_vector(PC::PolyhedralComplex)

Compute the vector $(f₀,f₁,f₂,...,f_{dim(PC))$` where $f_i$ is the number of
faces of $PC$ of dimension $i$.

# Examples
```jldoctest
julia> VR = [0 0; 1 0; -1 0; 0 1];

julia> IM = IncidenceMatrix([[1,2,4],[1,3,4]]);

julia> far_vertices = [2,3,4];

julia> PC = PolyhedralComplex(IM, VR, far_vertices);

julia> f_vector(PC)
3-element Vector{Int64}:
 1
 3
 2
```
"""
function f_vector(PC::PolyhedralComplex)
    ldim = lineality_dim(PC)
    f_vec=vcat(zeros(Int64, ldim), [length(polyhedra_of_dim(PC,i)) for i in ldim:dim(PC)])
    return f_vec
end


@doc Markdown.doc"""
    nrays(PC::PolyhedralComplex)

Return the number of rays of `PC`.

# Examples
```jldoctest
julia> VR = [0 0; 1 0; -1 0; 0 1];

julia> IM = IncidenceMatrix([[1,2,4],[1,3,4]]);

julia> far_vertices = [2,3,4];

julia> PC = PolyhedralComplex(IM, VR, far_vertices);

julia> nrays(PC)
3
```
"""
nrays(PC::PolyhedralComplex) = length(pm_object(PC).FAR_VERTICES)


@doc Markdown.doc"""
    nvertices(PC::PolyhedralComplex)

Return the number of vertices of `PC`.

# Examples
```jldoctest
julia> VR = [0 0; 1 0; -1 0; 0 1];

julia> IM = IncidenceMatrix([[1,2,4],[1,3,4]]);

julia> far_vertices = [2,3,4];

julia> PC = PolyhedralComplex(IM, VR, far_vertices);

julia> nvertices(PC)
1
```
"""
nvertices(PC::PolyhedralComplex) = pm_object(PC).N_VERTICES - nrays(PC)


@doc Markdown.doc"""
    npolyhedra(PC::PolyhedralComplex)

Return the total number of polyhedra in the polyhedral complex `PC`.

# Examples
```jldoctest
julia> VR = [0 0; 1 0; -1 0; 0 1];

julia> IM = IncidenceMatrix([[1,2,4],[1,3,4]]);

julia> far_vertices = [2,3,4];

julia> PC = PolyhedralComplex(IM, VR, far_vertices);

julia> npolyhedra(PC)
6
```
"""
npolyhedra(PC::PolyhedralComplex) = sum(f_vector(PC))

@doc Markdown.doc"""
    codim(PC::PolyhedralComplex)

Compute the codimension of a polyhedral complex.

# Examples
```
julia> VR = [0 0; 1 0; -1 0; 0 1];

julia> IM = IncidenceMatrix([[1,2],[1,3],[1,4]]);

julia> far_vertices = [2,3,4];

julia> PC = PolyhedralComplex(IM, VR, far_vertices)
A polyhedral complex in ambient dimension 2

julia> codim(PC)
1
```
"""
codim(PC::PolyhedralComplex) = ambient_dim(PC)-dim(PC)


@doc Markdown.doc"""
    isembedded(PC::PolyhedralComplex)

Returns true if `PC` is embedded, i.e. if its vertices can be computed as a
subset of some $\mathbb{R}^n$.

# Examples
```jldoctest
julia> VR = [0 0; 1 0; -1 0; 0 1];

julia> IM = IncidenceMatrix([[1,2],[1,3],[1,4]]);

julia> PC = PolyhedralComplex(IM, VR)
A polyhedral complex in ambient dimension 2

julia> isembedded(PC)
true
```
"""
function isembedded(PC::PolyhedralComplex)
    pmo = pm_object(PC)
    schedule = Polymake.call_method(pmo,:get_schedule,"VERTICES")
    return schedule != nothing
end
