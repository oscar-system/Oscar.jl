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
 Polytope in ambient dimension 3
 Polytope in ambient dimension 3
 Polytope in ambient dimension 3
 Polytope in ambient dimension 3
 Polytope in ambient dimension 3
 Polytope in ambient dimension 3
```
"""
function faces(P::Polyhedron{T}, face_dim::Int) where {T<:scalar_types}
  face_dim == dim(P) - 1 &&
    return SubObjectIterator{Polyhedron{T}}(P, _face_polyhedron_facet, n_facets(P))
  n = face_dim - length(lineality_space(P))
  n < 0 && return nothing
  pfaces = Polymake.to_one_based_indexing(Polymake.polytope.faces_of_dim(pm_object(P), n))
  farf = Polymake.to_one_based_indexing(pm_object(P).FAR_FACE)
  rfaces = Vector{Int}(filter(i -> !(pfaces[i] <= farf), range(1, length(pfaces); step=1)))
  return SubObjectIterator{Polyhedron{T}}(
    P, _face_polyhedron, length(rfaces), (f_dim=n, f_ind=rfaces)
  )
end

function _face_polyhedron(
  ::Type{Polyhedron{T}},
  P::Polyhedron{T},
  i::Base.Integer;
  f_dim::Int=-1,
  f_ind::Vector{Int64}=Vector{Int64}(),
) where {T<:scalar_types}
  pface = Polymake.to_one_based_indexing(
    Polymake.polytope.faces_of_dim(pm_object(P), f_dim)
  )[f_ind[i]]
  V = pm_object(P).VERTICES[collect(pface), :]
  L = pm_object(P).LINEALITY_SPACE
  PT = _scalar_type_to_polymake(T)
  return Polyhedron{T}(
    Polymake.polytope.Polytope{PT}(; VERTICES=V, LINEALITY_SPACE=L), coefficient_field(P)
  )
end

function _vertex_indices(
  ::Val{_face_polyhedron}, P::Polyhedron; f_dim=-1, f_ind::Vector{Int64}=Vector{Int64}()
)
  return IncidenceMatrix(
    collect.(
      Polymake.to_one_based_indexing(
        Polymake.polytope.faces_of_dim(pm_object(P), f_dim)[f_ind]
      )
    ),
  )[
    :, _vertex_indices(pm_object(P))
  ]
end

function _ray_indices(
  ::Val{_face_polyhedron}, P::Polyhedron; f_dim=-1, f_ind::Vector{Int64}=Vector{Int64}()
)
  return IncidenceMatrix(
    collect.(
      Polymake.to_one_based_indexing(
        Polymake.polytope.faces_of_dim(pm_object(P), f_dim)[f_ind]
      )
    ),
  )[
    :, _ray_indices(pm_object(P))
  ]
end

function _vertex_and_ray_indices(
  ::Val{_face_polyhedron}, P::Polyhedron; f_dim=-1, f_ind::Vector{Int64}=Vector{Int64}()
)
  return IncidenceMatrix(
    collect.(
      Polymake.to_one_based_indexing(
        Polymake.polytope.faces_of_dim(pm_object(P), f_dim)[f_ind]
      )
    ),
  )
end

_incidencematrix(::Val{_face_polyhedron}) = _vertex_and_ray_indices

function _face_polyhedron_facet(
  ::Type{Polyhedron{T}}, P::Polyhedron{T}, i::Base.Integer
) where {T<:scalar_types}
  pface = pm_object(P).VERTICES_IN_FACETS[_facet_index(pm_object(P), i), :]
  V = pm_object(P).VERTICES[collect(pface), :]
  L = pm_object(P).LINEALITY_SPACE
  PT = _scalar_type_to_polymake(T)
  return Polyhedron{T}(
    Polymake.polytope.Polytope{PT}(; VERTICES=V, LINEALITY_SPACE=L), coefficient_field(P)
  )
end

_vertex_indices(::Val{_face_polyhedron_facet}, P::Polyhedron) = vcat(
  pm_object(P).VERTICES_IN_FACETS[
    1:(_facet_at_infinity(pm_object(P)) - 1), _vertex_indices(pm_object(P))
  ],
  pm_object(P).VERTICES_IN_FACETS[
    (_facet_at_infinity(pm_object(P)) + 1):end, _vertex_indices(pm_object(P))
  ],
)

_ray_indices(::Val{_face_polyhedron_facet}, P::Polyhedron) = vcat(
  pm_object(P).VERTICES_IN_FACETS[
    1:(_facet_at_infinity(pm_object(P)) - 1), _ray_indices(pm_object(P))
  ],
  pm_object(P).VERTICES_IN_FACETS[
    (_facet_at_infinity(pm_object(P)) + 1):end, _ray_indices(pm_object(P))
  ],
)

_vertex_and_ray_indices(::Val{_face_polyhedron_facet}, P::Polyhedron) = vcat(
  pm_object(P).VERTICES_IN_FACETS[1:(_facet_at_infinity(pm_object(P)) - 1), :],
  pm_object(P).VERTICES_IN_FACETS[(_facet_at_infinity(pm_object(P)) + 1):end, :],
)

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

See also [`vertices`](@ref vertices(as::Type{PointVector{T}}, P::Polyhedron{T}) where {T<:scalar_types}) and [`lineality_space`](@ref lineality_space(P::Polyhedron{T}) where {T<:scalar_types}).

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
minimal_faces(P::Polyhedron{T}) where {T<:scalar_types} = minimal_faces(
  NamedTuple{
    (:base_points, :lineality_basis),
    Tuple{SubObjectIterator{PointVector{T}},SubObjectIterator{RayVector{T}}},
  },
  P,
)
function minimal_faces(
  ::Type{
    NamedTuple{
      (:base_points, :lineality_basis),
      Tuple{SubObjectIterator{PointVector{T}},SubObjectIterator{RayVector{T}}},
    },
  },
  P::Polyhedron{T},
) where {T<:scalar_types}
  return (base_points=_vertices(PointVector{T}, P), lineality_basis=lineality_space(P))
end
minimal_faces(::Type{<:PointVector}, P::Polyhedron) = _vertices(P)

@doc raw"""
    rays_modulo_lineality(as, P::Polyhedron)

Return the rays of the recession cone of `P` up to lineality as a `NamedTuple`
with two iterators. If `P` has lineality `L`, then the iterator
`rays_modulo_lineality` iterates over representatives of the rays of `P/L`.
The iterator `lineality_basis` gives a basis of the lineality space `L`.

See also [`rays`](@ref rays(as::Type{RayVector{T}}, P::Polyhedron{T}) where {T<:scalar_types}) and [`lineality_space`](@ref lineality_space(P::Polyhedron{T}) where {T<:scalar_types}).

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
rays_modulo_lineality(P::Polyhedron{T}) where {T<:scalar_types} = rays_modulo_lineality(
  NamedTuple{
    (:rays_modulo_lineality, :lineality_basis),
    Tuple{SubObjectIterator{RayVector{T}},SubObjectIterator{RayVector{T}}},
  },
  P,
)
function rays_modulo_lineality(
  ::Type{
    NamedTuple{
      (:rays_modulo_lineality, :lineality_basis),
      Tuple{SubObjectIterator{RayVector{T}},SubObjectIterator{RayVector{T}}},
    },
  },
  P::Polyhedron{T},
) where {T<:scalar_types}
  return (rays_modulo_lineality=_rays(P), lineality_basis=lineality_space(P))
end
rays_modulo_lineality(::Type{<:RayVector}, P::Polyhedron) = _rays(P)

@doc raw"""
    vertices([as::Type{T} = PointVector,] P::Polyhedron)

Return an iterator over the vertices of `P` in the format defined by `as`. The
vertices are defined to be the zero-dimensional faces, so if `P` has lineality,
there are no vertices, only minimal faces.

See also [`minimal_faces`](@ref minimal_faces(P::Polyhedron{T}) where {T<:scalar_types}) and [`rays`](@ref rays(as::Type{RayVector{T}}, P::Polyhedron{T}) where {T<:scalar_types}).

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

julia> point_matrix(vertices(P))
[-1   -1]
[ 2   -1]
[ 2    1]
[-1    2]
[ 1    2]
```
A half-space (here in $3$-space) has no vertices:
```jldoctest
julia> UH = polyhedron([-1 0 0], [0])
Polyhedron in ambient dimension 3

julia> vertices(UH)
0-element SubObjectIterator{PointVector{QQFieldElem}}
```
"""
vertices(as::Type{PointVector{T}}, P::Polyhedron{T}) where {T<:scalar_types} =
  lineality_dim(P) == 0 ? _vertices(as, P) : _empty_subobjectiterator(as, P)
_vertices(as::Type{PointVector{T}}, P::Polyhedron{T}) where {T<:scalar_types} =
  SubObjectIterator{as}(P, _vertex_polyhedron, length(_vertex_indices(pm_object(P))))

_vertex_polyhedron(
  U::Type{PointVector{T}}, P::Polyhedron{T}, i::Base.Integer
) where {T<:scalar_types} = point_vector(
  coefficient_field(P),
  @view pm_object(P).VERTICES[_vertex_indices(pm_object(P))[i], 2:end]
)::U

_point_matrix(::Val{_vertex_polyhedron}, P::Polyhedron; homogenized=false) =
  @view pm_object(P).VERTICES[_vertex_indices(pm_object(P)), (homogenized ? 1 : 2):end]

_matrix_for_polymake(::Val{_vertex_polyhedron}) = _point_matrix

vertices(::Type{<:PointVector}, P::Polyhedron{T}) where {T<:scalar_types} =
  vertices(PointVector{T}, P)
_vertices(::Type{<:PointVector}, P::Polyhedron{T}) where {T<:scalar_types} =
  _vertices(PointVector{T}, P)

_facet_indices(::Val{_vertex_polyhedron}, P::Polyhedron) =
  pm_object(P).FACETS_THRU_VERTICES[
    _vertex_indices(pm_object(P)), _facet_indices(pm_object(P))
  ]

_incidencematrix(::Val{_vertex_polyhedron}) = _facet_indices

function _facet_indices(P::Polymake.BigObject)
  vi = Polymake.get_attachment(P, "_facet_indices")
  if isnothing(vi)
    vi = Polymake.Vector{Polymake.to_cxx_type(Int64)}(
      [
        collect(1:(_facet_at_infinity(P) - 1))
        collect((_facet_at_infinity(P) + 1):size(P.FACETS, 1))
      ],
    )
    Polymake.attach(P, "_facet_indices", vi)
  end
  return vi
end

vertices(P::Polyhedron) = vertices(PointVector, P)
_vertices(P::Polyhedron) = _vertices(PointVector, P)

@doc raw"""
    n_rays(P::Polyhedron)

Return the number of rays of `P`, i.e. the number of rays of the recession cone
of `P`.

# Examples
The two-dimensional positive orthant has two rays.
```jldoctest
julia> PO = convex_hull([0 0],[1 0; 0 1])
Polyhedron in ambient dimension 2

julia> n_rays(PO)
2
```
The upper half-plane has no ray, since it has lineality.
```jldoctest
julia> UH = convex_hull([0 0],[0 1],[1 0]);

julia> n_rays(UH)
0
```
"""
n_rays(P::Polyhedron)::Int = lineality_dim(P) == 0 ? _n_rays(P) : 0
_n_rays(P::Polyhedron) = length(pm_object(P).FAR_FACE)

@doc raw"""
    n_vertices(P::Polyhedron)

Return the number of vertices of `P`.

# Examples
The 3-cube's number of vertices can be obtained with this input:
```jldoctest
julia> C = cube(3);

julia> n_vertices(C)
8
```
"""
n_vertices(P::Polyhedron)::Int = lineality_dim(P) == 0 ? _n_vertices(P) : 0
_n_vertices(P::Polyhedron) = size(pm_object(P).VERTICES, 1)::Int - _n_rays(P)

@doc raw"""
    rays([as::Type{T} = RayVector,] P::Polyhedron)

Return a minimal set of generators of the cone of unbounded directions of `P`
(i.e. its rays) in the format defined by `as`. The rays are defined to be the
one-dimensional faces of the recession cone, so if `P` has lineality, there
are no rays.

See also [`rays_modulo_lineality`](@ref rays_modulo_lineality(P::Polyhedron{T}) where {T<:scalar_types}), [`recession_cone`](@ref recession_cone(P::Polyhedron{T}) where {T<:scalar_types}) and [`vertices`](@ref vertices(as::Type{PointVector{T}}, P::Polyhedron{T}) where {T<:scalar_types}).

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
A half-space has no rays:
```
julia> UH = polyhedron([-1 0 0], [0])
Polyhedron in ambient dimension 3

julia> rays(UH)
0-element SubObjectIterator{RayVector{QQFieldElem}}
```
"""
rays(as::Type{RayVector{T}}, P::Polyhedron{T}) where {T<:scalar_types} =
  lineality_dim(P) == 0 ? _rays(as, P) : _empty_subobjectiterator(as, P)
_rays(as::Type{RayVector{T}}, P::Polyhedron{T}) where {T<:scalar_types} =
  SubObjectIterator{as}(P, _ray_polyhedron, length(_ray_indices(pm_object(P))))

_ray_polyhedron(
  U::Type{RayVector{T}}, P::Polyhedron{T}, i::Base.Integer
) where {T<:scalar_types} = ray_vector(
  coefficient_field(P), @view pm_object(P).VERTICES[_ray_indices(pm_object(P))[i], 2:end]
)::U

_facet_indices(::Val{_ray_polyhedron}, P::Polyhedron) =
  pm_object(P).FACETS_THRU_RAYS[_ray_indices(pm_object(P)), _facet_indices(pm_object(P))]

_vector_matrix(::Val{_ray_polyhedron}, P::Polyhedron; homogenized=false) =
  @view pm_object(P).VERTICES[_ray_indices(pm_object(P)), (homogenized ? 1 : 2):end]

_matrix_for_polymake(::Val{_ray_polyhedron}) = _vector_matrix

_incidencematrix(::Val{_ray_polyhedron}) = _facet_indices

rays(::Type{<:RayVector}, P::Polyhedron{T}) where {T<:scalar_types} = rays(RayVector{T}, P)
_rays(::Type{<:RayVector}, P::Polyhedron{T}) where {T<:scalar_types} =
  _rays(RayVector{T}, P)

rays(P::Polyhedron) = rays(RayVector, P)
_rays(P::Polyhedron) = _rays(RayVector, P)

@doc raw"""
    n_facets(P::Polyhedron)

Return the number of facets of `P`.

# Examples
The number of facets of the 5-dimensional cross polytope can be retrieved via
the following line:
```jldoctest
julia> n_facets(cross_polytope(5))
32
```
"""
function n_facets(P::Polyhedron)
  n = size(pm_object(P).FACETS, 1)::Int
  return n - (_facet_at_infinity(pm_object(P)) != n + 1)
end

@doc raw"""
    facets(as::Type{T} = AffineHalfspace, P::Polyhedron)

Return the facets of `P` in the format defined by `as`.

The allowed values for `as` are
* `Halfspace` (or its subtype `AffineHalfspace`),
* `Hyperplane` (or its subtype `AffineHyperplane`),
* `Polyhedron`,
* `Pair`.

# Examples
We can retrieve the six facets of the 3-dimensional cube this way:
```jldoctest
julia> C = cube(3);

julia> facets(Polyhedron, C)
6-element SubObjectIterator{Polyhedron{QQFieldElem}}:
 Polytope in ambient dimension 3
 Polytope in ambient dimension 3
 Polytope in ambient dimension 3
 Polytope in ambient dimension 3
 Polytope in ambient dimension 3
 Polytope in ambient dimension 3

julia> facets(Halfspace, C)
6-element SubObjectIterator{AffineHalfspace{QQFieldElem}} over the halfspaces of R^3 described by:
-x_1 <= 1
x_1 <= 1
-x_2 <= 1
x_2 <= 1
-x_3 <= 1
x_3 <= 1
```
"""
facets(
  as::Type{T}, P::Polyhedron{S}
) where {
  S<:scalar_types,
  T<:Union{AffineHalfspace{S},AffineHyperplane{S},Pair{R,S} where R,Polyhedron{S}},
} =
  SubObjectIterator{as}(P, _facet_polyhedron, n_facets(P))

function _facet_polyhedron(
  U::Type{AffineHalfspace{S}}, P::Polyhedron{S}, i::Base.Integer
) where {S<:scalar_types}
  h = decompose_hdata(view(pm_object(P).FACETS, [_facet_index(pm_object(P), i)], :))
  return affine_halfspace(coefficient_field(P), h[1], h[2][])::U
end
function _facet_polyhedron(
  U::Type{Pair{R,S}}, P::Polyhedron{S}, i::Base.Integer
) where {R,S<:scalar_types}
  f = coefficient_field(P)
  h = decompose_hdata(view(pm_object(P).FACETS, [_facet_index(pm_object(P), i)], :))
  return U(f.(view(h[1], :, :)), f(h[2][])) # view_broadcast
end
function _facet_polyhedron(
  ::Type{Polyhedron{T}}, P::Polyhedron{T}, i::Base.Integer
) where {T<:scalar_types}
  return Polyhedron{T}(
    Polymake.polytope.facet(pm_object(P), _facet_index(pm_object(P), i) - 1),
    coefficient_field(P),
  )
end
function _facet_polyhedron(
  U::Type{AffineHyperplane{S}}, P::Polyhedron{S}, i::Base.Integer
) where {S<:scalar_types}
  h = decompose_hdata(view(pm_object(P).FACETS, [_facet_index(pm_object(P), i)], :))
  return affine_hyperplane(coefficient_field(P), h[1], h[2][])::U
end

_affine_inequality_matrix(::Val{_facet_polyhedron}, P::Polyhedron) =
  -_remove_facet_at_infinity(pm_object(P))

_affine_matrix_for_polymake(::Val{_facet_polyhedron}) = _affine_inequality_matrix

_vertex_indices(::Val{_facet_polyhedron}, P::Polyhedron) = vcat(
  pm_object(P).VERTICES_IN_FACETS[
    1:(_facet_at_infinity(pm_object(P)) - 1), _vertex_indices(pm_object(P))
  ],
  pm_object(P).VERTICES_IN_FACETS[
    (_facet_at_infinity(pm_object(P)) + 1):end, _vertex_indices(pm_object(P))
  ],
)

_ray_indices(::Val{_facet_polyhedron}, P::Polyhedron) = vcat(
  pm_object(P).VERTICES_IN_FACETS[
    1:(_facet_at_infinity(pm_object(P)) - 1), _ray_indices(pm_object(P))
  ],
  pm_object(P).VERTICES_IN_FACETS[
    (_facet_at_infinity(pm_object(P)) + 1):end, _ray_indices(pm_object(P))
  ],
)

_vertex_and_ray_indices(::Val{_facet_polyhedron}, P::Polyhedron) = vcat(
  pm_object(P).VERTICES_IN_FACETS[1:(_facet_at_infinity(pm_object(P)) - 1), :],
  pm_object(P).VERTICES_IN_FACETS[(_facet_at_infinity(pm_object(P)) + 1):end, :],
)

_incidencematrix(::Val{_facet_polyhedron}) = _vertex_and_ray_indices

facets(::Type{<:Pair}, P::Polyhedron{T}) where {T<:scalar_types} =
  facets(Pair{Matrix{T},T}, P)

facets(::Type{Polyhedron}, P::Polyhedron{T}) where {T<:scalar_types} =
  facets(Polyhedron{T}, P)

facets(P::Polyhedron{T}) where {T<:scalar_types} = facets(AffineHalfspace{T}, P)

facets(::Type{<:Halfspace}, P::Polyhedron{T}) where {T<:scalar_types} =
  facets(AffineHalfspace{T}, P)
facets(::Type{<:Hyperplane}, P::Polyhedron{T}) where {T<:scalar_types} =
  facets(AffineHyperplane{T}, P)

function _facet_index(P::Polymake.BigObject, i::Base.Integer)
  i < _facet_at_infinity(P) && return i
  return i + 1
end

function _facet_at_infinity(P::Polymake.BigObject)
  fai = Polymake.get_attachment(P, "_facet_at_infinity")
  m = size(P.FACETS, 1)
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

_remove_facet_at_infinity(P::Polymake.BigObject) = view(
  P.FACETS,
  [
    collect(1:(_facet_at_infinity(P) - 1))
    collect((_facet_at_infinity(P) + 1):size(P.FACETS, 1))
  ],
  :,
)

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
volume(P::Polyhedron{T}) where {T<:scalar_types} =
  coefficient_field(P)((pm_object(P)).VOLUME)

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
lattice_volume(P::Polyhedron{QQFieldElem})::ZZRingElem =
  _assert_lattice(P) && pm_object(P).LATTICE_VOLUME

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
normalized_volume(P::Polyhedron) =
  coefficient_field(P)(factorial(dim(P)) * (pm_object(P)).VOLUME)

@doc raw"""
    castelnuovo_excess(P::Polyhedron)

For an arbitrary lattice polytope, Hibi [Hib94](@cite) proved that the normalized volume is always at least as large as a certain lattice point count.
This function returns the difference between those two numbers.

# Examples
```jldoctest
julia> castelnuovo_excess(cube(4))
154
```
"""
function castelnuovo_excess(P::Polyhedron)
  _assert_lattice(P)
  d = dim(P)
  l = ZZ(pm_object(P).N_LATTICE_POINTS)::ZZRingElem
  c = ZZ(pm_object(P).N_INTERIOR_LATTICE_POINTS)::ZZRingElem
  b = l - c
  e = (d * c + (d - 1) * b - d^2 + 2)
  # the normalized volume is an integer as P is lattice
  return numerator(normalized_volume(P)) - e
end

@doc raw"""
    is_castelnuovo(P::Polyhedron)

For an arbitrary lattice polytope, Hibi [Hib94](@cite) proved that the normalized volume is always at least as large as a certain lattice point count.
This function returns true if both numbers agree.

# Examples
```jldoctest
julia> is_castelnuovo(cube(2))
true

julia> is_castelnuovo(cube(4))
false
```
"""
function is_castelnuovo(P::Polyhedron)
  return castelnuovo_excess(P) == 0
end

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
    lattice_points(P::Polyhedron)

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
function lattice_points(P::Polyhedron; check::Bool=true)
  @check is_bounded(P) "Polyhedron not bounded"
  return SubObjectIterator{PointVector{ZZRingElem}}(
    P, _lattice_point, size(pm_object(P).LATTICE_POINTS_GENERATORS[1], 1)
  )
end

_lattice_point(
  T::Type{PointVector{ZZRingElem}}, P::Polyhedron, i::Base.Integer
) = point_vector(ZZ, @view pm_object(P).LATTICE_POINTS_GENERATORS[1][i, 2:end])::T

_point_matrix(::Val{_lattice_point}, P::Polyhedron; homogenized=false) =
  @view pm_object(P).LATTICE_POINTS_GENERATORS[1][:, (homogenized ? 1 : 2):end]

_matrix_for_polymake(::Val{_lattice_point}) = _point_matrix

@doc raw"""
    interior_lattice_points(P::Polyhedron)

Return the integer points contained in the interior of the bounded polyhedron
`P`.

# Examples
```jldoctest
julia> c = cube(3)
Polytope in ambient dimension 3

julia> interior_lattice_points(c)
1-element SubObjectIterator{PointVector{ZZRingElem}}:
 [0, 0, 0]

julia> matrix(ZZ, interior_lattice_points(c))
[0   0   0]
```
"""
function interior_lattice_points(P::Polyhedron)
  @req is_bounded(P) "Polyhedron not bounded"
  return SubObjectIterator{PointVector{ZZRingElem}}(
    P, _interior_lattice_point, size(pm_object(P).INTERIOR_LATTICE_POINTS, 1)
  )
end

_interior_lattice_point(
  T::Type{PointVector{ZZRingElem}}, P::Polyhedron, i::Base.Integer
) = point_vector(ZZ, @view pm_object(P).INTERIOR_LATTICE_POINTS[i, 2:end])::T

_point_matrix(::Val{_interior_lattice_point}, P::Polyhedron; homogenized=false) =
  if homogenized
    pm_object(P).INTERIOR_LATTICE_POINTS
  else
    @view pm_object(P).INTERIOR_LATTICE_POINTS[:, 2:end]
  end

_matrix_for_polymake(::Val{_interior_lattice_point}) = _point_matrix

@doc raw"""
    boundary_lattice_points(P::Polyhedron)

Return the integer points contained in the boundary of the bounded polyhedron
`P`.

# Examples
```jldoctest
julia> c = polarize(cube(3))
Polytope in ambient dimension 3

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
function boundary_lattice_points(P::Polyhedron)
  @req is_bounded(P) "Polyhedron not bounded"
  return SubObjectIterator{PointVector{ZZRingElem}}(
    P, _boundary_lattice_point, size(pm_object(P).BOUNDARY_LATTICE_POINTS, 1)
  )
end

_boundary_lattice_point(
  T::Type{PointVector{ZZRingElem}}, P::Polyhedron, i::Base.Integer
) = point_vector(ZZ, @view pm_object(P).BOUNDARY_LATTICE_POINTS[i, 2:end])::T

_point_matrix(::Val{_boundary_lattice_point}, P::Polyhedron; homogenized=false) =
  if homogenized
    pm_object(P).BOUNDARY_LATTICE_POINTS
  else
    @view pm_object(P).BOUNDARY_LATTICE_POINTS[:, 2:end]
  end

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
codim(P::Polyhedron) = ambient_dim(P) - dim(P)

@doc raw"""
    facet_sizes(P::Polyhedron{T})

Number of vertices in each facet. 

# Examples
```jldoctest
julia> p = johnson_solid(4) 
Polytope in ambient dimension 3 with EmbeddedAbsSimpleNumFieldElem type coefficients

julia> facet_sizes(p)
10-element Vector{Int64}:
 8
 4
 3
 4
 4
 3
 4
 3
 3
 4
```
"""
function facet_sizes(P::Polyhedron{T}) where {T<:scalar_types}
  im = vertex_indices(facets(P))
  return [length(row(im, i)) for i in 1:nrows(im)]
end

@doc raw"""
    vertex_sizes(P::Polyhedron{T})

Number of incident facets for each vertex.

# Examples
```jldoctest
julia> vertex_sizes(bipyramid(simplex(2)))
5-element Vector{Int64}:
 4
 4
 4
 3
 3
```
"""
function vertex_sizes(P::Polyhedron{T}) where {T<:scalar_types}
  pm_object(P).LINEALITY_DIM > 0 && return Vector{Int}()
  res = Vector{Int}(pm_object(P).VERTEX_SIZES)

  vertices = Polymake.to_one_based_indexing(
    Polymake.polytope.bounded_vertices(pm_object(P))
  )
  return res[collect(vertices)]
end

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
lineality_space(P::Polyhedron{T}) where {T<:scalar_types} =
  SubObjectIterator{RayVector{T}}(P, _lineality_polyhedron, lineality_dim(P))

_lineality_polyhedron(
  U::Type{RayVector{T}}, P::Polyhedron{T}, i::Base.Integer
) where {T<:scalar_types} =
  ray_vector(coefficient_field(P), @view pm_object(P).LINEALITY_SPACE[i, 2:end])::U

_generator_matrix(::Val{_lineality_polyhedron}, P::Polyhedron; homogenized=false) =
  homogenized ? pm_object(P).LINEALITY_SPACE : @view pm_object(P).LINEALITY_SPACE[:, 2:end]

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
2-element SubObjectIterator{AffineHyperplane{QQFieldElem}} over the hyperplanes of R^4 described by:
x_3 = 2
x_4 = 5
```
"""
affine_hull(P::Polyhedron{T}) where {T<:scalar_types} =
  SubObjectIterator{AffineHyperplane{T}}(P, _affine_hull, size(pm_object(P).AFFINE_HULL, 1))

function _affine_hull(
  U::Type{AffineHyperplane{T}}, P::Polyhedron{T}, i::Base.Integer
) where {T<:scalar_types}
  h = decompose_hdata(-view(pm_object(P).AFFINE_HULL, [i], :))
  return affine_hyperplane(coefficient_field(P), h[1], h[2][])::U
end

_affine_equation_matrix(::Val{_affine_hull}, P::Polyhedron) = pm_object(P).AFFINE_HULL

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
recession_cone(P::Polyhedron{T}) where {T<:scalar_types} =
  Cone{T}(Polymake.polytope.recession_cone(pm_object(P)), coefficient_field(P))

@doc raw"""
    ehrhart_polynomial(P::Polyhedron{QQFieldElem})

Compute the Ehrhart polynomial of `P`.

# Examples
```jldoctest
julia> c = cube(3)
Polytope in ambient dimension 3

julia> ehrhart_polynomial(c)
8*x^3 + 12*x^2 + 6*x + 1
```
"""
function ehrhart_polynomial(P::Polyhedron{QQFieldElem})
  R, x = polynomial_ring(QQ, :x; cached=false)
  return ehrhart_polynomial(R, P)
end

@doc raw"""
    ehrhart_polynomial(R::QQMPolyRing, P::Polyhedron{QQFieldElem})

Compute the Ehrhart polynomial of `P` and return it as a polynomial in `R`.

# Examples
```jldoctest
julia> R, x = polynomial_ring(QQ, :x)
(Univariate polynomial ring in x over QQ, x)

julia> c = cube(3)
Polytope in ambient dimension 3

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
Polytope in ambient dimension 3

julia> h_star_polynomial(c)
x^3 + 23*x^2 + 23*x + 1
```
"""
function h_star_polynomial(P::Polyhedron{QQFieldElem})
  R, x = polynomial_ring(QQ, :x; cached=false)
  return h_star_polynomial(R, P)
end

@doc raw"""
    h_star_polynomial(R::QQMPolyRing, P::Polyhedron)

Compute the $h^*$ polynomial of `P` and return it as a polynomial in `R`.

# Examples
```jldoctest
julia> R, x = polynomial_ring(QQ, :x)
(Univariate polynomial ring in x over QQ, x)

julia> c = cube(3)
Polytope in ambient dimension 3

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
Polytope in ambient dimension 3

julia> is_lattice_polytope(c)
true

julia> c = cube(3, 0, 4//3)
Polytope in ambient dimension 3

julia> is_lattice_polytope(c)
false
```
"""
is_lattice_polytope(P::Polyhedron{QQFieldElem}) =
  (is_bounded(P) && pm_object(P).LATTICE)::Bool

_assert_lattice(P::Polyhedron{QQFieldElem}) =
  is_lattice_polytope(P) ||
  throw(ArgumentError("This is only defined for lattice polytopes."))

@doc raw"""
    is_very_ample(P::Polyhedron{QQFieldElem})

Check whether `P` is very ample.

# Examples
```jldoctest
julia> c = cube(3)
Polytope in ambient dimension 3

julia> is_very_ample(c)
true

julia> P = convex_hull([0 0 0; 1 1 0; 1 0 1; 0 1 1])
Polyhedron in ambient dimension 3

julia> is_very_ample(P)
false
```
"""
is_very_ample(P::Polyhedron{QQFieldElem}) =
  _assert_lattice(P) && pm_object(P).VERY_AMPLE::Bool

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
Polytope in ambient dimension 3

julia> Q = cube(3,-1,2)
Polytope in ambient dimension 3

julia> issubset(P, Q)
true

julia> issubset(Q, P)
false
```
"""
Base.issubset(P::Polyhedron{T}, Q::Polyhedron{T}) where {T<:scalar_types} =
  Polymake.polytope.included_polyhedra(pm_object(P), pm_object(Q))::Bool

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
Base.in(v::AbstractVector, P::Polyhedron) =
  Polymake.polytope.contains(pm_object(P), coefficient_field(P).([1; v]))::Bool

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
Polytope in ambient dimension 3

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

Check whether `P` is simple, i.e., each vertex figure is a simplex.

# Examples
```jldoctest
julia> c = cube(2,0,1)
Polytope in ambient dimension 2

julia> is_simple(c)
true
```
"""
is_simple(P::Polyhedron) = pm_object(P).SIMPLE::Bool

@doc raw"""
    is_simplicial(P::Polyhedron)

Check whether `P` is simplicial, i.e., each proper face is a simplex.
"""
is_simplicial(P::Polyhedron) = pm_object(P).SIMPLICIAL::Bool

@doc raw"""
    is_neighborly(P::Polyhedron)

Check whether `P` is neighborly, i.e., if the dimension is $d$, each $\lfloor d/2 \rfloor$-subset of the vertices forms a face.
Neighborly polytopes in even dimension are necessarily simplicial.

# Examples

A 4-polytope is neighborly if and only if the vertex-edge graph is complete.

```jldoctest
julia> is_neighborly(cyclic_polytope(4,8))
true
```
"""
is_neighborly(P::Polyhedron) = pm_object(P).NEIGHBORLY::Bool

@doc raw"""
    is_cubical(P::Polyhedron)

Check whether `P` is cubical, i.e., each proper face is combinatorially equivalent to a cube.

# Examples

For details concerning the following construction see [JZ00](@cite).

```jldoctest
julia> Q = cube(2,-1,1); Q2 = cube(2,-2,2); P = convex_hull(product(Q,Q2), product(Q2,Q))
Polyhedron in ambient dimension 4

julia> is_cubical(P)
true
```
"""
is_cubical(P::Polyhedron) = pm_object(P).CUBICAL::Bool

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
    is_johnson_solid(P::Polyhedron)

Check whether `P` is a Johnson solid, i.e., a $3$-dimensional polytope with regular faces that is not vertex transitive.

See also [`johnson_solid`](@ref johnson_solid(index::Int)).

!!! note
    This will only recognize algebraically precise solids, i.e. no solids with approximate coordinates.

# Examples
```
julia> J = johnson_solid(37)
Polytope in ambient dimension 3 with EmbeddedAbsSimpleNumFieldElem type coefficients

julia> is_johnson_solid(J)
true
```
"""
is_johnson_solid(P::Polyhedron) = _is_3d_pol_reg_facets(P) && !is_vertex_transitive(P)

@doc raw"""
    is_archimedean_solid(P::Polyhedron)

Check whether `P` is an Archimedean solid, i.e., a $3$-dimensional vertex
transitive polytope with regular facets, but not a prism or antiprism.

See also [`archimedean_solid`](@ref archimedean_solid(s::String)).

!!! note
    This will only recognize algebraically precise solids, i.e. no solids with approximate coordinates.

# Examples
```jldoctest
julia> TO = archimedean_solid("truncated_octahedron")
Polytope in ambient dimension 3

julia> is_archimedean_solid(TO)
true

julia> T = tetrahedron()
Polytope in ambient dimension 3

julia> is_archimedean_solid(T)
false
```
"""
is_archimedean_solid(P::Polyhedron) =
  _is_3d_pol_reg_facets(P) && !_has_equal_facets(P) && is_vertex_transitive(P) &&
  !_is_prismic_or_antiprismic(P)

@doc raw"""
    is_platonic_solid(P::Polyhedron)

Check whether `P` is a Platonic solid.

See also [`platonic_solid`](@ref platonic_solid(s::String)).

!!! note
    This will only recognize algebraically precise solids, i.e. no solids with approximate coordinates.

# Examples
```jldoctest
julia> is_platonic_solid(cube(3))
true
```
"""
is_platonic_solid(P::Polyhedron) =
  _is_3d_pol_reg_facets(P) && _has_equal_facets(P) && is_vertex_transitive(P)

function _is_3d_pol_reg_facets(P::Polyhedron)
  # dimension
  dim(P) == 3 || return false

  # constant edge length
  pedges = faces(P, 1)

  edgelength = let v = vertices(pedges[1])
    _squared_distance(v[1], v[2])
  end

  for edge in pedges[2:end]
    v = vertices(edge)
    _squared_distance(v[1], v[2]) == edgelength || return false
  end

  pfacets = facets(Polyhedron, P)

  # facet vertices on circle
  for facet in pfacets
    n = n_vertices(facet)
    if n >= 4
      fverts = vertices(facet)
      m = sum(fverts)//n
      dist = _squared_distance(fverts[1], m)

      for v in fverts[2:end]
        _squared_distance(v, m) == dist || return false
      end
    end
  end
  return true
end

function _squared_distance(p::PointVector, q::PointVector)
  d = p - q
  return dot(d, d)
end

function _has_equal_facets(P::Polyhedron)
  return allequal(facet_sizes(P))
end

@doc raw"""
    is_vertex_transitive(P::Polyhedron)

Check whether `P` is vertex transitive.

# Examples
```jldoctest
julia> is_vertex_transitive(cube(3))
true

julia> is_vertex_transitive(pyramid(cube(2)))
false
```
"""
is_vertex_transitive(P::Polyhedron) =
  is_transitive(automorphism_group(P; action=:on_vertices))

function _is_prismic_or_antiprismic(P::Polyhedron)
  nvf = facet_sizes(P)
  # nvfs contains vector entries [i, j] where j is the amount of i-gonal facets of P
  nvfs = [[i, sum(nvf .== i)] for i in unique(nvf)]
  # an (anti-)prism can only have 2 different types of facets
  # n-gons and either squares/triangles (for prisms/antiprisms respectively)
  # if n == 1 this function will not be called either way
  length(nvfs) == 2 || return false
  # there are exactly 2 n-gons, and only n-gons can occur exactly twice
  a = findfirst(x -> x[2] == 2, nvfs)
  isnothing(a) && return false
  # with length(nvfs) == 2, b is the other index than a
  b = 3 - a
  m = nvfs[b][1]
  # the facets that are not n-gons have to be squares or triangles
  m == 3 || m == 4 || return false
  n = nvfs[a][1]
  # P has to have exactly n vertices and
  # the amount of squares needs to be n or the amount of triangles needs to be 2n
  2n == n_vertices(P) && (5 - m) * n == nvfs[b][2] || return false
  dg = dual_graph(P)
  ngon_is = findall(==(n), nvf)
  has_edge(dg, ngon_is...) && return false
  rem_vertex!(dg, ngon_is[2])
  rem_vertex!(dg, ngon_is[1])
  degs = [length(neighbors(dg, i)) for i in 1:n_vertices(dg)]
  degs == fill(2, n_vertices(dg)) && is_connected(dg) || return false
  dg = dual_graph(P)
  if m == 3
    bigdg = Polymake.graph.Graph(; ADJACENCY=dg.pm_graph)
    return bigdg.BIPARTITE
  else
    return true
  end
end

@doc raw"""
    f_vector(P::Polyhedron)

Return the vector $(f₀,f₁,f₂,...,f_{dim(P) - 1})$ where $f_i$ is the number of
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
function f_vector(P::Polyhedron)
  # the following differs from polymake's count in the unbounded case;
  # polymake takes the far face into account, too
  ldim = lineality_dim(P)
  f_vec = vcat(zeros(Int64, ldim), [length(faces(P, i)) for i in ldim:(dim(P) - 1)])
  return Vector{ZZRingElem}(f_vec)
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
function h_vector(P::Polyhedron)
  @req is_bounded(P) "defined for bounded polytopes only"
  return Vector{ZZRingElem}(pm_object(P).H_VECTOR)
end

@doc raw"""
    g_vector(P::Polyhedron)

Return the (toric) $g$-vector of a polytope.
Defined by $g_0 = 1$ and $g_k = h_k - h_{k-1}$, for $1 \leq k \leq \lceil (d+1)/2\rceil$ where $h$ is the $h$-vector and $d=\dim(P)$.
Undefined for unbounded polyhedra.

# Examples
```jldoctest
julia> g_vector(cross_polytope(3))
2-element Vector{ZZRingElem}:
 1
 2
```
"""
function g_vector(P::Polyhedron)
  @req is_bounded(P) "defined for bounded polytopes only"
  return Vector{ZZRingElem}(pm_object(P).G_VECTOR)
end

@doc raw"""
    relative_interior_point(P::Polyhedron)

Compute a point in the relative interior point of `P`, i.e. a point in `P` not
contained in any facet.

# Examples
The square $[-1,1]^3$ has the origin as a relative interior point.
```jldoctest
julia> square = cube(2)
Polytope in ambient dimension 2

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
relative_interior_point(P::Polyhedron{T}) where {T<:scalar_types} = point_vector(
  coefficient_field(P),
  view(dehomogenize(Polymake.common.dense(pm_object(P).REL_INT_POINT)), :),
)::PointVector{T} #view_broadcast

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
function support_function(P::Polyhedron{T}; convention=:max) where {T<:scalar_types}
  function h(ω::AbstractVector)
    lp = linear_program(P, ω; convention=convention)
    return solve_lp(lp)[1]
  end
  return h
end

_cmp_string(::Val{:lte}) = is_unicode_allowed() ? "≦" : "<="
_cmp_string(::Val{:eq}) = "="

@doc raw"""
    print_constraints([io = stdout,] A::AnyVecOrMat, b::AbstractVector; trivial = false, numbered = false, cmp = :lte)

Pretty print the constraints given by $P(A,b) = \{ x |  Ax ≤ b \}$.

# Optional & Keyword Arguments
- `io::IO`: Target `IO` where the  constraints are printed to.
- `trivial::Bool`: If `true`, include trivial inequalities.
- `numbered::Bool`: If `true`, the each constraint is printed with the index corresponding to the input `AnyVecOrMat`.
- `cmp::Symbol`: Defines the string used for the comparison sign; supports `:lte` (less than or equal) and `:eq` (equal).

Trivial inequalities are always counted for numbering, even when omitted.

# Examples
```jldoctest
julia> print_constraints([-1 0 4 5; 4 4 4 3; 1 0 0 0; 0 0 0 0; 0 0 0 0; 9 9 9 9], [0, 1, 2, 3, -4, 5]; numbered = true)
1: -x_1 + 4*x_3 + 5*x_4 <= 0
2: 4*x_1 + 4*x_2 + 4*x_3 + 3*x_4 <= 1
3: x_1 <= 2
5: 0 <= -4
6: 9*x_1 + 9*x_2 + 9*x_3 + 9*x_4 <= 5

julia> print_constraints([-1 0 4 5; 4 4 4 3; 1 0 0 0; 0 0 0 0; 0 0 0 0; 9 9 9 9], [0, 1, 2, 3, -4, 5]; trivial = true)
-x_1 + 4*x_3 + 5*x_4 <= 0
4*x_1 + 4*x_2 + 4*x_3 + 3*x_4 <= 1
x_1 <= 2
0 <= 3
0 <= -4
9*x_1 + 9*x_2 + 9*x_3 + 9*x_4 <= 5
```
"""
function print_constraints(
  io::IO,
  A::AnyVecOrMat,
  b::AbstractVector;
  trivial::Bool=false,
  numbered::Bool=false,
  cmp::Symbol=:lte,
)
  zero_char = is_unicode_allowed() ? '₀' : '0'
  us_ascii = is_unicode_allowed() ? "" : "_"
  for i in 1:length(b)
    terms = Vector{String}(undef, size(A)[2])
    first = true
    for j in 1:size(A)[2]
      if iszero(A[i, j])
        terms[j] = ""
      else
        if isone(A[i, j]) || isone(-A[i, j])
          terms[j] = if first
            string(
              isone(A[i, j]) ? "x" : "-x",
              us_ascii,
              reverse([zero_char + d for d in digits(j)])...,
            )
          else
            string(
              isone(A[i, j]) ? " + x" : " - x",
              us_ascii,
              reverse([zero_char + d for d in digits(j)])...,
            )
          end
        else
          terms[j] = if first
            string(
              _constraint_string(A[i, j]),
              "*x",
              us_ascii,
              reverse([zero_char + d for d in digits(j)])...,
            )
          else
            string(
              if A[i, j] < zero(A[i, j])
                string(" - ", _constraint_string(-A[i, j]))
              else
                string(" + ", _constraint_string(A[i, j]))
              end,
              "*x",
              us_ascii,
              reverse([zero_char + d for d in digits(j)])...,
            )
          end
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
    println(
      io,
      string(
        numbered ? string(i, ": ") : "", terms..., " ", _cmp_string(Val(cmp)), " ", b[i]
      ),
    )
  end
end

_constraint_string(x::Any) = string(x)
_constraint_string(x::QQFieldElem) = string(x)
_constraint_string(x::FieldElem) = string("(", x, ")")

@doc raw"""
    print_constraints([io = stdout,] P::Polyhedron; trivial = false, numbered = false)

Pretty print the constraints given by $P(A,b) = \{ x |  Ax ≤ b \}$.

# Optional & Keyword Arguments
- `io::IO`: Target `IO` where the  constraints are printed to.
- `trivial::Bool`: If `true`, include trivial inequalities.
- `numbered::Bool`: If `true`, the each constraint is printed with the index corresponding to the input `AnyVecOrMat`.

Trivial inequalities are always counted for numbering, even when omitted.

# Examples
The 3-cube is given by $-1 ≦ x_i ≦ 1 ∀ i ∈ \{1, 2, 3\}$.
```jldoctest
julia> print_constraints(cube(3))
-x_1 <= 1
x_1 <= 1
-x_2 <= 1
x_2 <= 1
-x_3 <= 1
x_3 <= 1
```
"""
print_constraints(io::IO, P::Polyhedron; trivial::Bool=false, numbered::Bool=false) =
  print_constraints(io, halfspace_matrix_pair(facets(P))...; trivial=trivial)

print_constraints(io::IO, H::Halfspace; trivial::Bool=false) =
  print_constraints(io, permutedims(normal_vector(H)), [negbias(H)]; trivial=trivial)

print_constraints(io::IO, H::Hyperplane; trivial::Bool=false) = print_constraints(
  io, permutedims(normal_vector(H)), [negbias(H)]; trivial=trivial, cmp=:eq
)

print_constraints(io::IO, H::SubObjectIterator{<:Halfspace}; numbered::Bool=false) =
  print_constraints(io, halfspace_matrix_pair(H)...; trivial=true, numbered=numbered)

print_constraints(io::IO, H::SubObjectIterator{<:Hyperplane}; numbered::Bool=false) =
  print_constraints(
    io, halfspace_matrix_pair(H)...; trivial=true, numbered=numbered, cmp=:eq
  )

# Default `io = stdout`
print_constraints(
  A::AnyVecOrMat,
  b::AbstractVector;
  trivial::Bool=false,
  numbered::Bool=false,
  cmp::Symbol=:lte,
) = print_constraints(stdout, A, b; trivial=trivial, numbered=numbered, cmp=cmp)

print_constraints(P::Polyhedron; trivial::Bool=false, numbered::Bool=false) =
  print_constraints(stdout, P; trivial=trivial, numbered=numbered)

print_constraints(H::Union{Halfspace,Hyperplane}; trivial::Bool=false) =
  print_constraints(stdout, H; trivial=trivial)

print_constraints(
  H::SubObjectIterator{<:Union{Halfspace,Hyperplane}}; numbered::Bool=false
) = print_constraints(stdout, H; numbered=numbered)

function Base.show(io::IO, H::Halfspace)
  n = length(normal_vector(H))
  if iszero(normal_vector(H)) && negbias(H) >= 0
    print(io, "The trivial half-space, R^$n")
  else
    print(io, "The half-space of R^$n described by\n")
    print_constraints(io, H)
  end
end

function Base.show(io::IO, H::Hyperplane)
  n = length(normal_vector(H))
  b = negbias(H)
  if iszero(b) && iszero(normal_vector(H))
    print(io, "The trivial hyperplane, R^$n")
  else
    print(io, "The hyperplane of R^$n described by\n")
    print_constraints(io, H)
  end
end

Base.show(io::IO, ::MIME"text/plain", H::SubObjectIterator{<:Union{Halfspace,Hyperplane}}) =
  show(io, H)

function Base.show(io::IO, H::SubObjectIterator{<:Halfspace})
  s = length(H)
  t = typeof(H)
  d = displaysize(io)[1] - 5
  print(io, "$s-element $t")
  if !isempty(H)
    n = length(normal_vector(H[1]))
    print(io, " over the halfspaces of R^$n described by:\n")
    if s < d
      print_constraints(io, H)
    else
      A, b = halfspace_matrix_pair(H)
      print_constraints(io, view(A, 1:floor(Int, d / 2), :), b[1:floor(Int, d / 2)])
      println(io, "⋮")
      print_constraints(
        io,
        A[(s - floor(Int, d / 2) + d % 2):end, :],
        b[(s - floor(Int, d / 2) + d % 2):end],
      )
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
    print(io, " over the hyperplanes of R^$n described by:\n")
    if s < d
      print_constraints(io, H)
    else
      A, b = halfspace_matrix_pair(H)
      print_constraints(io, A[1:floor(Int, d / 2), :], b[1:floor(Int, d / 2)]; cmp=:eq)
      println(io, "⋮")
      print_constraints(
        io,
        A[(s - floor(Int, d / 2) + d % 2):end, :],
        b[(s - floor(Int, d / 2) + d % 2):end];
        cmp=:eq,
      )
    end
  end
end
