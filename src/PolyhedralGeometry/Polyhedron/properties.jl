###############################################################################
###############################################################################
### Iterators
###############################################################################
###############################################################################

#TODO: take into account lineality space

@doc raw"""
    faces(P::Polyhedron, face_dim::Int)

Return an iterator over the faces of `P` of dimension `face_dim`.

# Examples
A `Vector` containing the six sides of the 3-dimensional cube can be obtained
via the following input:
```jldoctest
julia> F = faces(cube(3), 2)
6-element SubObjectIterator{Polyhedron{QQFieldElem}}:
 Polyhedron in ambient dimension 3
 Polyhedron in ambient dimension 3
 Polyhedron in ambient dimension 3
 Polyhedron in ambient dimension 3
 Polyhedron in ambient dimension 3
 Polyhedron in ambient dimension 3
```
"""
function faces(P::Polyhedron{T}, face_dim::Int) where T<:scalar_types
    face_dim == dim(P) - 1 && return SubObjectIterator{Polyhedron{T}}(pm_object(P), _face_polyhedron_facet, nfacets(P))
    n = face_dim - length(lineality_space(P))
    n < 0 && return nothing
    pfaces = Polymake.to_one_based_indexing(Polymake.polytope.faces_of_dim(pm_object(P), n))
    farf = Polymake.to_one_based_indexing(pm_object(P).FAR_FACE)
    rfaces = Vector{Int}(filter(i->!(pfaces[i] <= farf), range(1,length(pfaces); step=1)))
    return SubObjectIterator{Polyhedron{T}}(pm_object(P), _face_polyhedron, length(rfaces), (f_dim = n, f_ind = rfaces))
end

function _face_polyhedron(::Type{Polyhedron{T}}, P::Polymake.BigObject, i::Base.Integer; f_dim::Int = -1, f_ind::Vector{Int64} = Vector{Int64}()) where T<:scalar_types
    pface = Polymake.to_one_based_indexing(Polymake.polytope.faces_of_dim(P, f_dim))[f_ind[i]]
    return Polyhedron{T}(Polymake.polytope.Polytope{scalar_type_to_polymake[T]}(VERTICES = P.VERTICES[collect(pface),:], LINEALITY_SPACE = P.LINEALITY_SPACE))
end

function _vertex_indices(::Val{_face_polyhedron}, P::Polymake.BigObject; f_dim = -1, f_ind::Vector{Int64} = Vector{Int64}())
    return IncidenceMatrix(collect.(Polymake.to_one_based_indexing(Polymake.polytope.faces_of_dim(P, f_dim)[f_ind])))[:, _vertex_indices(P)]
end

function _ray_indices(::Val{_face_polyhedron}, P::Polymake.BigObject; f_dim = -1, f_ind::Vector{Int64} = Vector{Int64}())
    return IncidenceMatrix(collect.(Polymake.to_one_based_indexing(Polymake.polytope.faces_of_dim(P, f_dim)[f_ind])))[:, _ray_indices(P)]
end

function _vertex_and_ray_indices(::Val{_face_polyhedron}, P::Polymake.BigObject; f_dim = -1, f_ind::Vector{Int64} = Vector{Int64}())
    return IncidenceMatrix(collect.(Polymake.to_one_based_indexing(Polymake.polytope.faces_of_dim(P, f_dim)[f_ind])))
end

_incidencematrix(::Val{_face_polyhedron}) = _vertex_and_ray_indices

function _face_polyhedron_facet(::Type{Polyhedron{T}}, P::Polymake.BigObject, i::Base.Integer) where T<:scalar_types
    pface = P.VERTICES_IN_FACETS[_facet_index(P, i), :]
    return Polyhedron{T}(Polymake.polytope.Polytope{scalar_type_to_polymake[T]}(VERTICES = P.VERTICES[collect(pface),:], LINEALITY_SPACE = P.LINEALITY_SPACE))
end

_vertex_indices(::Val{_face_polyhedron_facet}, P::Polymake.BigObject) = vcat(P.VERTICES_IN_FACETS[1:(_facet_at_infinity(P) - 1), _vertex_indices(P)], P.VERTICES_IN_FACETS[(_facet_at_infinity(P) + 1):end, _vertex_indices(P)])

_ray_indices(::Val{_face_polyhedron_facet}, P::Polymake.BigObject) = vcat(P.VERTICES_IN_FACETS[1:(_facet_at_infinity(P) - 1), _ray_indices(P)], P.VERTICES_IN_FACETS[(_facet_at_infinity(P) + 1):end, _ray_indices(P)])

_vertex_and_ray_indices(::Val{_face_polyhedron_facet}, P::Polymake.BigObject) = vcat(P.VERTICES_IN_FACETS[1:(_facet_at_infinity(P) - 1), :], P.VERTICES_IN_FACETS[(_facet_at_infinity(P) + 1):end, :])

_incidencematrix(::Val{_face_polyhedron_facet}) = _vertex_and_ray_indices

function _isray(P::Polyhedron, i::Base.Integer)
    return in(i, _ray_indices(pm_object(P)))
end

function _vertex_indices(P::Polymake.BigObject)
    vi = Polymake.get_attachment(P, "_vertex_indices")
    if isnothing(vi)
        A = P.VERTICES
        vi = Polymake.Vector{Polymake.to_cxx_type(Int64)}(findall(!iszero, view(A, :, 1)))
        Polymake.attach(P, "_vertex_indices", vi)
    end
    return vi
end

_ray_indices(P::Polymake.BigObject) = collect(Polymake.to_one_based_indexing(P.FAR_FACE))

function _polymake_to_oscar_vertex_index(P::Polymake.BigObject, i::Base.Integer)
    return i - sum((>).(i, P.FAR_FACE))
end

function _polymake_to_oscar_vertex_index(P::Polymake.BigObject, v::AbstractVector)
    return [_polymake_to_oscar_vertex_index(P, v[i]) for i in 1:length(v)]
end

function _polymake_to_oscar_ray_index(P::Polymake.BigObject, i::Base.Integer)
    return sum((<).(i, P.FAR_FACE))
end

function _polymake_to_oscar_ray_index(P::Polymake.BigObject, v::AbstractVector)
    return [_polymake_to_oscar_ray_index(P, v[i]) for i in 1:length(v)]
end


@doc raw"""
    minimal_faces(as, P::Polyhedron)

Return the minimal faces of a polyhedron as a `NamedTuple` with two iterators.
For a polyhedron without lineality, the `base_points` are the vertices. If `P`
has lineality `L`, then every minimal face is an affine translation `p+L`,
where `p` is only unique modulo `L`. The return type is a dict, the key
`:base_points` gives an iterator over such `p`, and the key `:lineality_basis`
lets one access a basis for the lineality space `L` of `P`.

# Examples
The polyhedron `P` is just a line through the origin:
```jldoctest
julia> P = convex_hull([0 0], nothing, [1 0])
Polyhedron in ambient dimension 2

julia> lineality_dim(P)
1

julia> vertices(P)
0-element SubObjectIterator{PointVector{QQFieldElem}}

julia> minimal_faces(P)
(base_points = PointVector{QQFieldElem}[[0, 0]], lineality_basis = RayVector{QQFieldElem}[[1, 0]])
```
"""
minimal_faces(P::Polyhedron{T}) where T<:scalar_types = minimal_faces(NamedTuple{(:base_points, :lineality_basis), Tuple{SubObjectIterator{PointVector{T}}, SubObjectIterator{RayVector{T}}}}, P)
function minimal_faces(as::Type{NamedTuple{(:base_points, :lineality_basis), Tuple{SubObjectIterator{PointVector{T}}, SubObjectIterator{RayVector{T}}}}}, P::Polyhedron{T}) where T<:scalar_types
    return (
        base_points = _vertices(PointVector{T}, P),
        lineality_basis = lineality_space(P)
    )
end
minimal_faces(as::Type{PointVector{T}}, P::Polyhedron{T}) where T<:scalar_types = _vertices(PointVector{T}, P)



@doc raw"""
    rays_modulo_lineality(as, P::Polyhedron)

Return the rays of the recession cone of `P` up to lineality as a `NamedTuple`
with two iterators. If `P` has lineality `L`, then the iterator
`rays_modulo_lineality` iterates over representatives of the rays of `P/L`.
The iterator `lineality_basis` gives a basis of the lineality space `L`.

# Examples
```jldoctest
julia> P = convex_hull([0 0 1], [0 1 0], [1 0 0])
Polyhedron in ambient dimension 3

julia> rmlP = rays_modulo_lineality(P)
(rays_modulo_lineality = RayVector{QQFieldElem}[[0, 1, 0]], lineality_basis = RayVector{QQFieldElem}[[1, 0, 0]])

julia> rmlP.rays_modulo_lineality
1-element SubObjectIterator{RayVector{QQFieldElem}}:
 [0, 1, 0]

julia> rmlP.lineality_basis
1-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0, 0]
```
"""
rays_modulo_lineality(P::Polyhedron{T}) where T<:scalar_types = rays_modulo_lineality(NamedTuple{(:rays_modulo_lineality, :lineality_basis), Tuple{SubObjectIterator{RayVector{T}}, SubObjectIterator{RayVector{T}}}}, P)
function rays_modulo_lineality(as::Type{NamedTuple{(:rays_modulo_lineality, :lineality_basis), Tuple{SubObjectIterator{RayVector{T}}, SubObjectIterator{RayVector{T}}}}}, P::Polyhedron) where T<:scalar_types
    return (
        rays_modulo_lineality = _rays(P),
        lineality_basis = lineality_space(P)
    )
end
rays_modulo_lineality(as::Type{RayVector}, P::Polyhedron) = _rays(P)


@doc raw"""
    vertices(as, P)

Return an iterator over the vertices of `P` in the format defined by `as`. The
vertices are defined to be the zero-dimensional faces, so if `P` has lineality,
there are no vertices.

Optional arguments for `as` include
* `PointVector`.

# Examples
The following code computes the vertices of the Minkowski sum of a triangle and
a square:
```jldoctest
julia> P = simplex(2) + cube(2);

julia> vertices(PointVector, P)
5-element SubObjectIterator{PointVector{QQFieldElem}}:
 [-1, -1]
 [2, -1]
 [2, 1]
 [-1, 2]
 [1, 2]
```
"""
vertices(as::Type{PointVector{T}}, P::Polyhedron) where T<:scalar_types = lineality_dim(P) == 0 ? _vertices(as, P) : _empty_subobjectiterator(as, pm_object(P))
_vertices(as::Type{PointVector{T}}, P::Polyhedron) where T<:scalar_types = SubObjectIterator{as}(pm_object(P), _vertex_polyhedron, length(_vertex_indices(pm_object(P))))

_vertex_polyhedron(::Type{PointVector{T}}, P::Polymake.BigObject, i::Base.Integer) where T<:scalar_types = PointVector{T}(@view P.VERTICES[_vertex_indices(P)[i], 2:end])

_point_matrix(::Val{_vertex_polyhedron}, P::Polymake.BigObject; homogenized=false) = @view P.VERTICES[_vertex_indices(P), (homogenized ? 1 : 2):end]

_matrix_for_polymake(::Val{_vertex_polyhedron}) = _point_matrix

vertices(::Type{PointVector}, P::Polyhedron{T}) where T<:scalar_types = vertices(PointVector{T}, P)
_vertices(::Type{PointVector}, P::Polyhedron{T}) where T<:scalar_types = _vertices(PointVector{T}, P)

_facet_indices(::Val{_vertex_polyhedron}, P::Polymake.BigObject)=P.FACETS_THRU_VERTICES[_vertex_indices(P),_facet_indices(P)]

_incidencematrix(::Val{_vertex_polyhedron}) = _facet_indices

function _facet_indices(P::Polymake.BigObject)
    vi = Polymake.get_attachment(P, "_facet_indices")
    if isnothing(vi)
        vi = Polymake.Vector{Polymake.to_cxx_type(Int64)}([collect(1:(_facet_at_infinity(P) - 1)); collect((_facet_at_infinity(P) + 1):size(P.FACETS, 1))])
        Polymake.attach(P, "_facet_indices", vi);
    end
    return vi;
end

@doc raw"""
    vertices(P::Polyhedron)

Return an iterator over the vertices of a polyhedron `P` as points.

# Examples
The following code computes the vertices of the Minkowski sum of a triangle and
a square:
```jldoctest
julia> P = simplex(2) + cube(2);

julia> vertices(P)
5-element SubObjectIterator{PointVector{QQFieldElem}}:
 [-1, -1]
 [2, -1]
 [2, 1]
 [-1, 2]
 [1, 2]
```
"""
vertices(P::Polyhedron) = vertices(PointVector, P)
_vertices(P::Polyhedron) = _vertices(PointVector, P)


@doc raw"""
    nrays(P::Polyhedron)

Return the number of rays of `P`, i.e. the number of rays of the recession cone
of `P`.

# Examples
The two-dimensional positive orthant has two rays.
```jldoctest
julia> PO = convex_hull([0 0],[1 0; 0 1])
Polyhedron in ambient dimension 2

julia> nrays(PO)
2
```
The upper half-plane has no ray, since it has lineality.
```jldoctest
julia> UH = convex_hull([0 0],[0 1],[1 0]);

julia> nrays(UH)
0
```
"""
nrays(P::Polyhedron)::Int = lineality_dim(P) == 0 ? _nrays(P) : 0
_nrays(P::Polyhedron) = length(pm_object(P).FAR_FACE)

@doc raw"""
    nvertices(P::Polyhedron)

Return the number of vertices of `P`.

# Examples
The 3-cube's number of vertices can be obtained with this input:
```jldoctest
julia> C = cube(3);

julia> nvertices(C)
8
```
"""
nvertices(P::Polyhedron)::Int = lineality_dim(P) == 0 ? _nvertices(P) : 0
_nvertices(P::Polyhedron) = size(pm_object(P).VERTICES, 1)::Int - _nrays(P)


@doc raw"""
    rays(as::Type{T} = RayVector, P::Polyhedron)

Return a minimal set of generators of the cone of unbounded directions of `P`
(i.e. its rays) in the format defined by `as`.

Optional arguments for `as` include
* `RayVector`.

# Examples
We can verify that the positive orthant of the plane is generated by the two
rays in positive unit direction:
```jldoctest
julia> PO = convex_hull([0 0], [1 0; 0 1]);

julia> rays(RayVector, PO)
2-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0]
 [0, 1]
```
"""
rays(as::Type{RayVector{T}}, P::Polyhedron) where T<:scalar_types = lineality_dim(P) == 0 ? _rays(as, P) : _empty_subobjectiterator(as, pm_object(P))
_rays(as::Type{RayVector{T}}, P::Polyhedron) where T<:scalar_types = SubObjectIterator{as}(pm_object(P), _ray_polyhedron, length(_ray_indices(pm_object(P))))

_ray_polyhedron(::Type{RayVector{T}}, P::Polymake.BigObject, i::Base.Integer) where T<:scalar_types = RayVector{T}(@view P.VERTICES[_ray_indices(P)[i], 2:end])

_facet_indices(::Val{_ray_polyhedron}, P::Polymake.BigObject)=P.FACETS_THRU_RAYS[_ray_indices(P),_facet_indices(P)]

_vector_matrix(::Val{_ray_polyhedron}, P::Polymake.BigObject; homogenized=false) = @view P.VERTICES[_ray_indices(P), (homogenized ? 1 : 2):end]

_matrix_for_polymake(::Val{_ray_polyhedron}) = _vector_matrix

_incidencematrix(::Val{_ray_polyhedron}) = _facet_indices

rays(::Type{RayVector}, P::Polyhedron{T}) where T<:scalar_types = rays(RayVector{T}, P)
_rays(::Type{RayVector}, P::Polyhedron{T}) where T<:scalar_types = _rays(RayVector{T}, P)

@doc raw"""
    rays(P::Polyhedron)

Return minimal set of generators of the cone of unbounded directions of `P`
(i.e. its rays) as points.

# Examples
We can verify that the positive orthant of the plane is generated by the two
rays in positive unit direction:
```jldoctest
julia> PO = convex_hull([0 0], [1 0; 0 1]);

julia> rays(PO)
2-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0]
 [0, 1]

julia> matrix(QQ, rays(PO))
[1   0]
[0   1]

julia> matrix(ZZ, rays(PO))
[1   0]
[0   1]
```
"""
rays(P::Polyhedron) = rays(RayVector, P)
_rays(P::Polyhedron) = _rays(RayVector, P)

@doc raw"""
    nfacets(P::Polyhedron)

Return the number of facets of `P`.

# Examples
The number of facets of the 5-dimensional cross polytope can be retrieved via
the following line:
```jldoctest
julia> nfacets(cross_polytope(5))
32
```
"""
function nfacets(P::Polyhedron)
    n = size(pm_object(P).FACETS, 1)::Int
    return n - (_facet_at_infinity(pm_object(P)) != n + 1)
end

@doc raw"""
    facets(as::Type{T} = AffineHalfspace, P::Polyhedron)

Return the facets of `P` in the format defined by `as`.

The allowed values for `as` are
* `Halfspace`,
* `Polyhedron`,
* `Pair`.

# Examples
We can retrieve the six facets of the 3-dimensional cube this way:
```jldoctest
julia> C = cube(3);

julia> facets(Polyhedron, C)
6-element SubObjectIterator{Polyhedron{QQFieldElem}}:
 Polyhedron in ambient dimension 3
 Polyhedron in ambient dimension 3
 Polyhedron in ambient dimension 3
 Polyhedron in ambient dimension 3
 Polyhedron in ambient dimension 3
 Polyhedron in ambient dimension 3

julia> facets(Halfspace, C)
6-element SubObjectIterator{AffineHalfspace{QQFieldElem}} over the Halfspaces of R^3 described by:
-x₁ ≦ 1
x₁ ≦ 1
-x₂ ≦ 1
x₂ ≦ 1
-x₃ ≦ 1
x₃ ≦ 1
```
"""
facets(as::Type{T}, P::Polyhedron{S}) where {R, S<:scalar_types, T<:Union{AffineHalfspace{S}, Pair{R, S}, Polyhedron{S}}} = SubObjectIterator{as}(pm_object(P), _facet_polyhedron, nfacets(P))

function _facet_polyhedron(::Type{T}, P::Polymake.BigObject, i::Base.Integer) where {R, S<:scalar_types, T<:Union{AffineHalfspace{S}, Pair{R, S}}}
    h = decompose_hdata(view(P.FACETS, [_facet_index(P, i)], :))
    return T(h[1], h[2][])
end
function _facet_polyhedron(::Type{Polyhedron{T}}, P::Polymake.BigObject, i::Base.Integer) where T<:scalar_types
    h = decompose_hdata(view(P.FACETS, [_facet_index(P, i)], :))
    return polyhedron(T, h[1], h[2][])
end

_affine_inequality_matrix(::Val{_facet_polyhedron}, P::Polymake.BigObject) = -_remove_facet_at_infinity(P)

_affine_matrix_for_polymake(::Val{_facet_polyhedron}) = _affine_inequality_matrix

_vertex_indices(::Val{_facet_polyhedron}, P::Polymake.BigObject) = vcat(P.VERTICES_IN_FACETS[1:(_facet_at_infinity(P) - 1), _vertex_indices(P)], P.VERTICES_IN_FACETS[(_facet_at_infinity(P) + 1):end, _vertex_indices(P)])

_ray_indices(::Val{_facet_polyhedron}, P::Polymake.BigObject) = vcat(P.VERTICES_IN_FACETS[1:(_facet_at_infinity(P) - 1), _ray_indices(P)], P.VERTICES_IN_FACETS[(_facet_at_infinity(P) + 1):end, _ray_indices(P)])

_vertex_and_ray_indices(::Val{_facet_polyhedron}, P::Polymake.BigObject) = vcat(P.VERTICES_IN_FACETS[1:(_facet_at_infinity(P) - 1), :], P.VERTICES_IN_FACETS[(_facet_at_infinity(P) + 1):end, :])

_incidencematrix(::Val{_facet_polyhedron}) = _vertex_and_ray_indices

facets(::Type{Pair}, P::Polyhedron{T}) where T<:scalar_types = facets(Pair{Matrix{T}, T}, P)

facets(::Type{Polyhedron}, P::Polyhedron{T}) where T<:scalar_types = facets(Polyhedron{T}, P)

@doc raw"""
    facets(P::Polyhedron)

Return the facets of `P` as halfspaces.

# Examples
We can retrieve the six facets of the 3-dimensional cube this way:
```jldoctest
julia> C = cube(3);

julia> facets(C)
6-element SubObjectIterator{AffineHalfspace{QQFieldElem}} over the Halfspaces of R^3 described by:
-x₁ ≦ 1
x₁ ≦ 1
-x₂ ≦ 1
x₂ ≦ 1
-x₃ ≦ 1
x₃ ≦ 1
```
"""
facets(P::Polyhedron{T}) where T<:scalar_types = facets(AffineHalfspace{T}, P)

facets(::Type{Halfspace}, P::Polyhedron{T}) where T<:scalar_types = facets(AffineHalfspace{T}, P)

facets(::Type{AffineHalfspace}, P::Polyhedron{T}) where T<:scalar_types = facets(AffineHalfspace{T}, P)

function _facet_index(P::Polymake.BigObject, i::Base.Integer)
    i < _facet_at_infinity(P) && return i
    return i + 1
end

function _facet_at_infinity(P::Polymake.BigObject)
    fai = Polymake.get_attachment(P, "_facet_at_infinity")
    m = size(P.FACETS,1)
    if isnothing(fai)
        i = 1
        while i <= m
            _is_facet_at_infinity(view(P.FACETS, i, :)) && break
            i += 1
        end
        fai = i
        Polymake.attach(P, "_facet_at_infinity", fai)
    end
    return fai::Int64
end

_is_facet_at_infinity(v::AbstractVector) = v[1] >= 0 && iszero(v[2:end])

_remove_facet_at_infinity(P::Polymake.BigObject) = view(P.FACETS, [collect(1:(_facet_at_infinity(P) - 1)); collect((_facet_at_infinity(P) + 1):size(P.FACETS, 1))], :)

###############################################################################
###############################################################################
### Access properties
###############################################################################
###############################################################################

###############################################################################
## Scalar properties
###############################################################################
@doc raw"""
    lineality_dim(P::Polyhedron)

Return the dimension of the lineality space, i.e. the dimension of the largest
affine subspace contained in `P`.

# Examples
Polyhedron with one lineality direction.
```jldoctest
julia> C = convex_hull([0 0], [1 0], [1 1])
Polyhedron in ambient dimension 2

julia> lineality_dim(C)
1
```
"""
lineality_dim(P::Polyhedron) = pm_object(P).LINEALITY_DIM::Int


@doc raw"""
    volume(P::Polyhedron)

Return the (Euclidean) volume of `P`.

# Examples
```jldoctest
julia> C = cube(2);

julia> volume(C)
4
```
"""
volume(P::Polyhedron{T}) where T<:scalar_types = convert(T, (pm_object(P)).VOLUME)

volume(P::Polyhedron{nf_elem}) = convert(nf_scalar, pm_object(P).VOLUME)


@doc raw"""
    lattice_volume(P::Polyhedron{QQFieldElem})

Return the lattice volume of `P`.

# Examples
```jldoctest
julia> C = cube(2);

julia> lattice_volume(C)
8
```
"""
lattice_volume(P::Polyhedron{QQFieldElem})::ZZRingElem = _assert_lattice(P) && pm_object(P).LATTICE_VOLUME


@doc raw"""
    normalized_volume(P::Polyhedron)

Return the (normalized) volume of `P`.

# Examples
```jldoctest
julia> C = cube(2);

julia> normalized_volume(C)
8
```
"""
normalized_volume(P::Polyhedron{T}) where T<:scalar_types = convert(T, factorial(dim(P))*(pm_object(P)).VOLUME)

normalized_volume(P::Polyhedron{nf_elem}) = convert(nf_scalar, factorial(dim(P))*(pm_object(P)).VOLUME)


@doc raw"""
    dim(P::Polyhedron)

Return the dimension of `P`.

# Examples
```jldoctest
julia> V = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];

julia> P = convex_hull(V);

julia> dim(P)
2
```
"""
dim(P::Polyhedron) = Polymake.polytope.dim(pm_object(P))::Int


@doc raw"""
    lattice_points(P::Polyhedron{QQFieldElem})

Return the integer points contained in the bounded polyhedron `P`.

# Examples
```jldoctest
julia> S = 2 * simplex(2);

julia> lattice_points(S)
6-element SubObjectIterator{PointVector{ZZRingElem}}:
 [0, 0]
 [0, 1]
 [0, 2]
 [1, 0]
 [1, 1]
 [2, 0]

julia> matrix(ZZ, lattice_points(S))
[0   0]
[0   1]
[0   2]
[1   0]
[1   1]
[2   0]
```
"""
function lattice_points(P::Polyhedron{QQFieldElem})
    @req pm_object(P).BOUNDED "Polyhedron not bounded"
    return SubObjectIterator{PointVector{ZZRingElem}}(pm_object(P), _lattice_point, size(pm_object(P).LATTICE_POINTS_GENERATORS[1], 1))
end

_lattice_point(::Type{PointVector{ZZRingElem}}, P::Polymake.BigObject, i::Base.Integer) = PointVector{ZZRingElem}(@view P.LATTICE_POINTS_GENERATORS[1][i, 2:end])

_point_matrix(::Val{_lattice_point}, P::Polymake.BigObject; homogenized=false) = @view P.LATTICE_POINTS_GENERATORS[1][:, (homogenized ? 1 : 2):end]

_matrix_for_polymake(::Val{_lattice_point}) = _point_matrix


@doc raw"""
    interior_lattice_points(P::Polyhedron{QQFieldElem})

Return the integer points contained in the interior of the bounded polyhedron
`P`.

# Examples
```jldoctest
julia> c = cube(3)
Polyhedron in ambient dimension 3

julia> interior_lattice_points(c)
1-element SubObjectIterator{PointVector{ZZRingElem}}:
 [0, 0, 0]

julia> matrix(ZZ, interior_lattice_points(c))
[0   0   0]
```
"""
function interior_lattice_points(P::Polyhedron{QQFieldElem})
    @req pm_object(P).BOUNDED "Polyhedron not bounded"
    return SubObjectIterator{PointVector{ZZRingElem}}(pm_object(P), _interior_lattice_point, size(pm_object(P).INTERIOR_LATTICE_POINTS, 1))
end

_interior_lattice_point(::Type{PointVector{ZZRingElem}}, P::Polymake.BigObject, i::Base.Integer) = PointVector{ZZRingElem}(@view P.INTERIOR_LATTICE_POINTS[i, 2:end])

_point_matrix(::Val{_interior_lattice_point}, P::Polymake.BigObject; homogenized=false) = homogenized ? P.INTERIOR_LATTICE_POINTS : @view P.INTERIOR_LATTICE_POINTS[:, 2:end]

_matrix_for_polymake(::Val{_interior_lattice_point}) = _point_matrix

@doc raw"""
    boundary_lattice_points(P::Polyhedron{QQFieldElem})

Return the integer points contained in the boundary of the bounded polyhedron
`P`.

# Examples
```jldoctest
julia> c = polarize(cube(3))
Polyhedron in ambient dimension 3

julia> boundary_lattice_points(c)
6-element SubObjectIterator{PointVector{ZZRingElem}}:
 [-1, 0, 0]
 [0, -1, 0]
 [0, 0, -1]
 [0, 0, 1]
 [0, 1, 0]
 [1, 0, 0]

julia> matrix(ZZ, boundary_lattice_points(c))
[-1    0    0]
[ 0   -1    0]
[ 0    0   -1]
[ 0    0    1]
[ 0    1    0]
[ 1    0    0]
```
"""
function boundary_lattice_points(P::Polyhedron{QQFieldElem})
    @req pm_object(P).BOUNDED "Polyhedron not bounded"
    return SubObjectIterator{PointVector{ZZRingElem}}(pm_object(P), _boundary_lattice_point, size(pm_object(P).BOUNDARY_LATTICE_POINTS, 1))
end

_boundary_lattice_point(::Type{PointVector{ZZRingElem}}, P::Polymake.BigObject, i::Base.Integer) = PointVector{ZZRingElem}(@view P.BOUNDARY_LATTICE_POINTS[i, 2:end])

_point_matrix(::Val{_boundary_lattice_point}, P::Polymake.BigObject; homogenized=false) = homogenized ? P.BOUNDARY_LATTICE_POINTS : @view P.BOUNDARY_LATTICE_POINTS[:, 2:end]

_matrix_for_polymake(::Val{_boundary_lattice_point}) = _point_matrix

@doc raw"""
    ambient_dim(P::Polyhedron)

Return the ambient dimension of `P`.

# Examples
```jldoctest
julia> V = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];

julia> P = convex_hull(V);

julia> ambient_dim(P)
3
```
"""
ambient_dim(P::Polyhedron) = Polymake.polytope.ambient_dim(pm_object(P))::Int


@doc raw"""
    codim(P::Polyhedron)

Return the codimension of `P`.

# Examples
```jldoctest
julia> V = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];

julia> P = convex_hull(V);

julia> codim(P)
1
```
"""
codim(P::Polyhedron) = ambient_dim(P)-dim(P)

###############################################################################
## Points properties
###############################################################################


# Previously: This implementation is not correct. Ask Taylor.
# Taylor: lineality space generators always look like [0, v] so
#  v is a natural output.
@doc raw"""
    lineality_space(P::Polyhedron)

Return a matrix whose row span is the lineality space of `P`.

# Examples
Despite not being reflected in this construction of the upper half-plane,
its lineality in $x$-direction is recognized:
```jldoctest
julia> UH = convex_hull([0 0],[0 1; 1 0; -1 0]);

julia> lineality_space(UH)
1-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0]
```
"""
lineality_space(P::Polyhedron{T}) where T<:scalar_types = SubObjectIterator{RayVector{T}}(pm_object(P), _lineality_polyhedron, lineality_dim(P))

_lineality_polyhedron(::Type{RayVector{T}}, P::Polymake.BigObject, i::Base.Integer) where T<:scalar_types = RayVector{T}(@view P.LINEALITY_SPACE[i, 2:end])

_generator_matrix(::Val{_lineality_polyhedron}, P::Polymake.BigObject; homogenized=false) = homogenized ? P.LINEALITY_SPACE : @view P.LINEALITY_SPACE[:, 2:end]

_matrix_for_polymake(::Val{_lineality_polyhedron}) = _generator_matrix


@doc raw"""
    affine_hull(P::Polytope)

Return the (affine) hyperplanes generating the affine hull of `P`.

# Examples
This triangle in $\mathbb{R}^4$ is contained in the plane defined by
$P = \{ (x_1, x_2, x_3, x_4) | x_3 = 2 ∧ x_4 = 5 \}$.
```jldoctest
julia> t = convex_hull([0 0 2 5; 1 0 2 5; 0 1 2 5]);

julia> affine_hull(t)
2-element SubObjectIterator{AffineHyperplane{QQFieldElem}} over the Hyperplanes of R^4 described by:
x₃ = 2
x₄ = 5
```
"""
affine_hull(P::Polyhedron{T}) where T<:scalar_types = SubObjectIterator{AffineHyperplane{T}}(pm_object(P), _affine_hull, size(pm_object(P).AFFINE_HULL, 1))

function _affine_hull(::Type{AffineHyperplane{T}}, P::Polymake.BigObject, i::Base.Integer) where T<:scalar_types
    h = decompose_hdata(-view(P.AFFINE_HULL, [i], :))
    return AffineHyperplane{T}(h[1], h[2][])
end

_affine_equation_matrix(::Val{_affine_hull}, P::Polymake.BigObject) = P.AFFINE_HULL

_affine_matrix_for_polymake(::Val{_affine_hull}) = _affine_equation_matrix


@doc raw"""
    recession_cone(P::Polyhedron)

Return the recession cone of `P`.

# Examples
```jldoctest
julia> P = polyhedron([1 -2; -1 1; -1 0; 0 -1],[2,1,1,1]);

julia> vertices(P)
3-element SubObjectIterator{PointVector{QQFieldElem}}:
 [0, -1]
 [-1, 0]
 [-1, -1]

julia> recession_cone(P)
Polyhedral cone in ambient dimension 2

julia> rays(recession_cone(P))
2-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 1//2]
 [1, 1]
```
"""
recession_cone(P::Polyhedron{T}) where T<:scalar_types = Cone{T}(Polymake.polytope.recession_cone(pm_object(P)))


@doc raw"""
    ehrhart_polynomial(P::Polyhedron{QQFieldElem})

Compute the Ehrhart polynomial of `P`.

# Examples
```jldoctest
julia> c = cube(3)
Polyhedron in ambient dimension 3

julia> ehrhart_polynomial(c)
8*x^3 + 12*x^2 + 6*x + 1
```
"""
function ehrhart_polynomial(P::Polyhedron{QQFieldElem})
    R, x = polynomial_ring(QQ, "x")
    return ehrhart_polynomial(R, P)
end


@doc raw"""
    ehrhart_polynomial(R::QQMPolyRing, P::Polyhedron{QQFieldElem})

Compute the Ehrhart polynomial of `P` and return it as a polynomial in `R`.

# Examples
```jldoctest
julia> R, x = polynomial_ring(QQ, "x")
(Univariate polynomial ring in x over QQ, x)

julia> c = cube(3)
Polyhedron in ambient dimension 3

julia> ehrhart_polynomial(R, c)
8*x^3 + 12*x^2 + 6*x + 1
```
"""
function ehrhart_polynomial(R::QQPolyRing, P::Polyhedron{QQFieldElem})
  _assert_lattice(P)
  coeffs = Polymake.polytope.ehrhart_polynomial_coeff(pm_object(P))
  return (R)(Vector{QQFieldElem}(coeffs))
end


@doc raw"""
    h_star_polynomial(P::Polyhedron)

Compute the $h^*$ polynomial of `P`.

# Examples
```jldoctest
julia> c = cube(3)
Polyhedron in ambient dimension 3

julia> h_star_polynomial(c)
x^3 + 23*x^2 + 23*x + 1
```
"""
function h_star_polynomial(P::Polyhedron{QQFieldElem})
    R, x = polynomial_ring(QQ, "x")
    return h_star_polynomial(R, P)
end


@doc raw"""
    h_star_polynomial(R::QQMPolyRing, P::Polyhedron)

Compute the $h^*$ polynomial of `P` and return it as a polynomial in `R`.

# Examples
```jldoctest
julia> R, x = polynomial_ring(QQ, "x")
(Univariate polynomial ring in x over QQ, x)

julia> c = cube(3)
Polyhedron in ambient dimension 3

julia> h_star_polynomial(R, c)
x^3 + 23*x^2 + 23*x + 1
```
"""
function h_star_polynomial(R::QQPolyRing, P::Polyhedron{QQFieldElem})
    coeffs = pm_object(P).H_STAR_VECTOR
    return (R)(Vector{QQFieldElem}(coeffs))
end

###############################################################################
## Boolean properties
###############################################################################
@doc raw"""
    is_lattice_polytope(P::Polyhedron{QQFieldElem})

Check whether `P` is a lattice polytope, i.e. it is bounded and has integral vertices.

# Examples
```jldoctest
julia> c = cube(3)
Polyhedron in ambient dimension 3

julia> is_lattice_polytope(c)
true

julia> c = cube(3, 0, 4//3)
Polyhedron in ambient dimension 3

julia> is_lattice_polytope(c)
false
```
"""
is_lattice_polytope(P::Polyhedron{QQFieldElem}) = (is_bounded(P) && pm_object(P).LATTICE)::Bool

_assert_lattice(P::Polyhedron{QQFieldElem}) = is_lattice_polytope(P) ||
  throw(ArgumentError("This is only defined for lattice polytopes."))


@doc raw"""
    is_very_ample(P::Polyhedron{QQFieldElem})

Check whether `P` is very ample.

# Examples
```jldoctest
julia> c = cube(3)
Polyhedron in ambient dimension 3

julia> is_very_ample(c)
true

julia> P = convex_hull([0 0 0; 1 1 0; 1 0 1; 0 1 1])
Polyhedron in ambient dimension 3

julia> is_very_ample(P)
false
```
"""
is_very_ample(P::Polyhedron{QQFieldElem}) = _assert_lattice(P) && pm_object(P).VERY_AMPLE::Bool


@doc raw"""
    is_feasible(P::Polyhedron)

Check whether `P` is feasible, i.e. non-empty.

# Examples
```jldoctest
julia> P = polyhedron([1 -1; -1 1; -1 0; 0 -1],[-1,-1,1,1]);

julia> is_feasible(P)
false
```
"""
is_feasible(P::Polyhedron) = pm_object(P).FEASIBLE::Bool


@doc raw"""
    issubset(P::Polyhedron, Q::Polyhedron)

Check whether `P` is a subset of the polyhedron `Q`.

# Examples
```jldoctest
julia> P = cube(3,0,1)
Polyhedron in ambient dimension 3

julia> Q = cube(3,-1,2)
Polyhedron in ambient dimension 3

julia> issubset(P, Q)
true

julia> issubset(Q, P)
false
```
"""
Base.issubset(P::Polyhedron{T}, Q::Polyhedron{T}) where T<:scalar_types = Polymake.polytope.included_polyhedra(pm_object(P), pm_object(Q))::Bool


@doc raw"""
    in(v::AbstractVector, P::Polyhedron)

Check whether the vector `v` is contained in the polyhedron `P`.

# Examples
The positive orthant only contains vectors with non-negative entries:
```jldoctest
julia> PO = polyhedron([-1 0; 0 -1], [0, 0]);

julia> [1, 2] in PO
true

julia> [1, -2] in PO
false
```
"""
Base.in(v::AbstractVector, P::Polyhedron) = Polymake.polytope.contains(pm_object(P), [1; v])::Bool


@doc raw"""
    is_smooth(P::Polyhedron{QQFieldElem})

Check whether `P` is smooth.

# Examples
A cube is always smooth.
```jldoctest
julia> C = cube(8);

julia> is_smooth(C)
true
```
"""
is_smooth(P::Polyhedron{QQFieldElem}) = _assert_lattice(P) && pm_object(P).SMOOTH::Bool


@doc raw"""
    is_normal(P::Polyhedron{QQFieldElem})

Check whether `P` is normal.

# Examples
The 3-cube is normal.
```jldoctest
julia> C = cube(3)
Polyhedron in ambient dimension 3

julia> is_normal(C)
true
```
But this pyramid is not:
```jldoctest
julia> P = convex_hull([0 0 0; 0 1 1; 1 1 0; 1 0 1]);

julia> is_normal(P)
false
```
"""
is_normal(P::Polyhedron{QQFieldElem}) = _assert_lattice(P) && pm_object(P).NORMAL::Bool


@doc raw"""
    is_bounded(P::Polyhedron)

Check whether `P` is bounded.

# Examples
```jldoctest
julia> P = polyhedron([1 -3; -1 1; -1 0; 0 -1],[1,1,1,1]);

julia> is_bounded(P)
false
```
"""
is_bounded(P::Polyhedron) = pm_object(P).BOUNDED::Bool


@doc raw"""
    is_simple(P::Polyhedron)

Check whether `P` is simple.

# Examples
```jldoctest
julia> c = cube(2,0,1)
Polyhedron in ambient dimension 2

julia> is_simple(c)
true
```
"""
is_simple(P::Polyhedron) = pm_object(P).SIMPLE::Bool


@doc raw"""
    is_simplicial(P::Polyhedron)

Check whether `P` is simplicial.
"""
is_simplicial(P::Polyhedron) = pm_object(P).SIMPLICIAL::Bool


@doc raw"""
    is_fulldimensional(P::Polyhedron)

Check whether `P` is full-dimensional.

# Examples
```jldoctest
julia> V = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];

julia> is_fulldimensional(convex_hull(V))
false
```
"""
is_fulldimensional(P::Polyhedron) = pm_object(P).FULL_DIM::Bool


@doc raw"""
    f_vector(P::Polyhedron)

Return the vector $(f₀,f₁,f₂,...,f_{(dim(P)-1))$` where $f_i$ is the number of
faces of $P$ of dimension $i$.

# Examples
Here we compute the f-vector of the 5-cube:
```jldoctest
julia> f_vector(cube(5))
5-element Vector{ZZRingElem}:
 32
 80
 80
 40
 10
```
"""
function f_vector(P::Polyhedron)::Vector{ZZRingElem}
    # the following differs from polymake's count in the unbounded case;
    # polymake takes the far face into account, too
    ldim = lineality_dim(P)
    f_vec=vcat(zeros(Int64, ldim), [length(faces(P,i)) for i in ldim:dim(P)-1])
    return f_vec
end

@doc raw"""
    h_vector(P::Polyhedron)

Return the (toric) h-vector of a polytope.
For simplicial polytopes this is a linear transformation of the f-vector.
Undefined for unbounded polyhedra.

# Examples
```jldoctest
julia> h_vector(cross_polytope(3))
4-element Vector{ZZRingElem}:
 1
 3
 3
 1
```
"""
function h_vector(P::Polyhedron)::Vector{ZZRingElem}
    @req is_bounded(P) "defined for bounded polytopes only"
    return pm_object(P).H_VECTOR
end


@doc raw"""
    g_vector(P::Polyhedron)

Return the (toric) $g$-vector of a polytope.
Defined by $g_0 = 1 $ and $g_k = h_k - h_{k-1}$, for $1 \leq k \leq \lceil (d+1)/2\rceil$ where $h$ is the $h$-vector and $d=\dim(P)$.
Undefined for unbounded polyhedra.

# Examples
```jldoctest
julia> g_vector(cross_polytope(3))
2-element Vector{ZZRingElem}:
 1
 2
```
"""
function g_vector(P::Polyhedron)::Vector{ZZRingElem}
    @req is_bounded(P) "defined for bounded polytopes only"
    return pm_object(P).G_VECTOR
end


@doc raw"""
    relative_interior_point(P::Polyhedron)

Compute a point in the relative interior point of `P`, i.e. a point in `P` not
contained in any facet.

# Examples
The square $[-1,1]^3$ has the origin as a relative interior point.
```jldoctest
julia> square = cube(2)
Polyhedron in ambient dimension 2

julia> relative_interior_point(square)
2-element PointVector{QQFieldElem}:
 0
 0

julia> vertices(square)
4-element SubObjectIterator{PointVector{QQFieldElem}}:
 [-1, -1]
 [1, -1]
 [-1, 1]
 [1, 1]

julia> matrix(QQ, vertices(square))
[-1   -1]
[ 1   -1]
[-1    1]
[ 1    1]
```
"""
relative_interior_point(P::Polyhedron{T}) where T<:scalar_types = PointVector{T}(dehomogenize(Polymake.common.dense(pm_object(P).REL_INT_POINT)))


@doc raw"""
    support_function(P::Polyhedron; convention::Symbol = :max)

Produce a function $h(ω) = max\{dot(x,ω)\ |\ x \in P\}$. $max$ may be changed
to $min$ by setting `convention = :min`.

# Examples
```jldoctest
julia> P = cube(3) + simplex(3);

julia> φ = support_function(P);

julia> φ([1,2,3])
9

julia> ψ = support_function(P, convention = :min);

julia> ψ([1,2,3])
-6
```
"""
function support_function(P::Polyhedron{T}; convention = :max) where T<:scalar_types
    function h(ω::AbstractVector)
        lp=linear_program(P,ω; convention = convention)
        return solve_lp(lp)[1]
    end
    return h
end


@doc raw"""
    print_constraints(A::AnyVecOrMat, b::AbstractVector; trivial::Bool = false, numbered::Bool = false)

Pretty print the constraints given by $P(A,b) = \{ x |  Ax ≤ b \}$.

Trivial inequalities are counted but omitted. They are included if `trivial` is
set to `true`.

# Examples
```jldoctest
julia> print_constraints([-1 0 4 5; 4 4 4 3; 1 0 0 0; 0 0 0 0; 0 0 0 0; 9 9 9 9], [0, 1, 2, 3, -4, 5]; numbered = true)
1: -x₁ + 4*x₃ + 5*x₄ ≦ 0
2: 4*x₁ + 4*x₂ + 4*x₃ + 3*x₄ ≦ 1
3: x₁ ≦ 2
5: 0 ≦ -4
6: 9*x₁ + 9*x₂ + 9*x₃ + 9*x₄ ≦ 5

julia> print_constraints([-1 0 4 5; 4 4 4 3; 1 0 0 0; 0 0 0 0; 0 0 0 0; 9 9 9 9], [0, 1, 2, 3, -4, 5]; trivial = true)
-x₁ + 4*x₃ + 5*x₄ ≦ 0
4*x₁ + 4*x₂ + 4*x₃ + 3*x₄ ≦ 1
x₁ ≦ 2
0 ≦ 3
0 ≦ -4
9*x₁ + 9*x₂ + 9*x₃ + 9*x₄ ≦ 5
```
"""
function print_constraints(A::AnyVecOrMat, b::AbstractVector; trivial::Bool = false, numbered::Bool = false, io::IO = stdout, cmp::String = "≦")
    for i in 1:length(b)
        terms = Vector{String}(undef, size(A)[2])
        first = true
        for j in 1:size(A)[2]
            if iszero(A[i, j])
                terms[j] = ""
            else
                if isone(A[i, j]) || isone(-A[i, j])
                    terms[j] = first ? string(isone(A[i, j]) ? "x" : "-x" , ['₀'+ d for d in digits(j)]...) :
                        string(isone(A[i, j]) ? " + x" : " - x", ['₀'+ d for d in digits(j)]...)
                else
                    terms[j] = first ? string(_constraint_string(A[i, j]), "*x", ['₀'+ d for d in digits(j)]...) :
                        string(A[i, j] < 0 ? string(" - ", _constraint_string(-A[i, j])) : string(" + ", _constraint_string(A[i, j])), "*x", ['₀'+ d for d in digits(j)]...)
                end
                first = false
            end
        end
        if first
            if b[i] >= 0 && !trivial
                continue
            end
            terms[1] = "0"
        end
        println(io, string(numbered ? string(i, ": ") : "", terms..., " ", cmp, " ", b[i]))
    end
end

_constraint_string(x::Any) = string(x)

_constraint_string(x::nf_elem) = string("(", x, ")")

@doc raw"""
    print_constraints(P::Polyhedron; trivial::Bool = false, numbered::Bool = false)

Pretty print the constraints given by $P(A,b) = \{ x |  Ax ≤ b \}$.

Trivial inequalities are counted but omitted. They are included if `trivial` is
set to `true`.

# Examples
The 3-cube is given by $-1 ≦ x_i ≦ 1 ∀ i ∈ \{1, 2, 3\}$.
```jldoctest
julia> print_constraints(cube(3))
-x₁ ≦ 1
x₁ ≦ 1
-x₂ ≦ 1
x₂ ≦ 1
-x₃ ≦ 1
x₃ ≦ 1
```
"""
print_constraints(P::Polyhedron; trivial::Bool = false, numbered::Bool = false, io::IO = stdout) = print_constraints(halfspace_matrix_pair(facets(P))...; trivial = trivial, io = io)

print_constraints(H::Halfspace; trivial::Bool = false, io::IO = stdout) = print_constraints(hcat(normal_vector(H)...), [negbias(H)]; trivial = trivial, io = io)

print_constraints(H::Hyperplane; trivial::Bool = false, io::IO = stdout) = print_constraints(hcat(normal_vector(H)...), [negbias(H)]; trivial = trivial, io = io, cmp = "=")

print_constraints(H::SubObjectIterator{<:Halfspace}; numbered::Bool = false, io::IO = stdout) = print_constraints(halfspace_matrix_pair(H)...; trivial = true, numbered = numbered, io = io)

print_constraints(H::SubObjectIterator{<:Hyperplane}; numbered::Bool = false, io::IO = stdout) = print_constraints(halfspace_matrix_pair(H)...; trivial = true, numbered = numbered, io = io, cmp = "=")

function Base.show(io::IO, H::Halfspace)
    n = length(normal_vector(H))
    if iszero(normal_vector(H)) && negbias(H) >= 0
        print(io, "The trivial half-space, R^$n")
    else
        print(io, "The half-space of R^$n described by\n")
        print_constraints(H; io=io)
    end
end

function Base.show(io::IO, H::Hyperplane)
    n = length(normal_vector(H))
    b = negbias(H)
    if iszero(b) && iszero(normal_vector(H))
        print(io, "The trivial hyperplane, R^$n")
    else
        print(io, "The hyperplane of R^$n described by\n")
        print_constraints(H; io = io)
    end
end

Base.show(io::IO, ::MIME"text/plain", H::SubObjectIterator{<:Union{Halfspace, Hyperplane}}) = show(io, H)

function Base.show(io::IO, H::SubObjectIterator{<:Halfspace})
    s = length(H)
    t = typeof(H)
    d = displaysize(io)[1] - 5
    print(io, "$s-element $t")
    if !isempty(H)
        n = length(normal_vector(H[1]))
        print(io, " over the Halfspaces of R^$n described by:\n")
        if s < d
            print_constraints(H; io = io)
        else
            A, b = halfspace_matrix_pair(H)
            print_constraints(view(A, 1:floor(Int, d/2), :), b[1:floor(Int, d/2)]; io = io)
            println(io, "⋮")
            print_constraints(A[(s - floor(Int, d/2) + d%2):end, :], b[(s - floor(Int, d/2) + d%2):end]; io = io)
        end
    end
end

function Base.show(io::IO, H::SubObjectIterator{<:Hyperplane})
    s = length(H)
    t = typeof(H)
    d = displaysize(io)[1] - 5
    print(io, "$s-element $t")
    if !isempty(H)
        n = length(normal_vector(H[1]))
        print(io, " over the Hyperplanes of R^$n described by:\n")
        if s < d
            print_constraints(H; io = io)
        else
            A, b = halfspace_matrix_pair(H)
            print_constraints(A[1:floor(Int, d/2), :], b[1:floor(Int, d/2)]; io = io, cmp = "=")
            println(io, "⋮")
            print_constraints(A[(s - floor(Int, d/2) + d%2):end, :], b[(s - floor(Int, d/2) + d%2):end]; io = io, cmp = "=")
        end
    end
end
