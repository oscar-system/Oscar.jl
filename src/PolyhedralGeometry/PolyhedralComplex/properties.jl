@doc raw"""
    ambient_dim(PC::PolyhedralComplex)

Return the ambient dimension of `PC`.

# Examples
```jldoctest
julia> IM = incidence_matrix([[1,2,3],[1,3,4]])
2×4 IncidenceMatrix
[1, 2, 3]
[1, 3, 4]


julia> V = [0 0; 1 0; 1 1; 0 1]
4×2 Matrix{Int64}:
 0  0
 1  0
 1  1
 0  1

julia> PC = polyhedral_complex(IM, V)
Polyhedral complex in ambient dimension 2

julia> ambient_dim(PC)
2
```
"""
ambient_dim(PC::PolyhedralComplex) = Polymake.fan.ambient_dim(pm_object(PC))::Int

@doc raw"""
    vertices([as::Type{T} = PointVector,] PC::PolyhedralComplex)

Return an iterator over the vertices of `PC` in the format defined by `as`. The
vertices are defined to be the zero-dimensional faces, so if `P` has lineality,
there are no vertices, only minimal faces.

See also [`minimal_faces`](@ref minimal_faces(PC::PolyhedralComplex{T}) where {T<:scalar_types}) and [`rays`](@ref rays(PC::PolyhedralComplex{T}) where {T<:scalar_types}).

Optional arguments for `as` include
* `PointVector`.

# Examples
```jldoctest
julia> IM = incidence_matrix([[1,2,3],[1,3,4]]);

julia> V = [0 0; 1 0; 1 1; 0 1];

julia> PC = polyhedral_complex(IM, V)
Polyhedral complex in ambient dimension 2

julia> vertices(PC)
4-element SubObjectIterator{PointVector{QQFieldElem}}:
 [0, 0]
 [1, 0]
 [1, 1]
 [0, 1]

julia> matrix(QQ, vertices(PointVector, PC))
[0   0]
[1   0]
[1   1]
[0   1]
```
The following complex has no vertices:
```
julia> IM = incidence_matrix([[1,2,3],[1,3,4]]);

julia> VR = [0 0 0; 1 0 0; 0 1 0; -1 0 0];

julia> far_vertices = [2,3,4];

julia> L = [0 0 1];

julia> PC = polyhedral_complex(IM, VR, far_vertices, L)
Polyhedral complex in ambient dimension 3

julia> vertices(PC)
0-element SubObjectIterator{PointVector{QQFieldElem}}
```
"""
vertices(as::Type{PointVector{T}}, PC::PolyhedralComplex{T}) where {T<:scalar_types} =
  lineality_dim(PC) == 0 ? _vertices(as, PC) : _empty_subobjectiterator(as, PC)
_vertices(as::Type{PointVector{T}}, PC::PolyhedralComplex{T}) where {T<:scalar_types} =
  SubObjectIterator{as}(PC, _vertex_complex, length(_vertex_indices(pm_object(PC))))

vertices(as::Type{<:PointVector}, PC::PolyhedralComplex{T}) where {T<:scalar_types} =
  vertices(PointVector{T}, PC)

_vertex_complex(
  U::Type{PointVector{T}}, PC::PolyhedralComplex{T}, i::Base.Integer
) where {T<:scalar_types} = point_vector(
  coefficient_field(PC),
  @view pm_object(PC).VERTICES[_vertex_indices(pm_object(PC))[i], 2:end]
)::U

_point_matrix(::Val{_vertex_complex}, PC::PolyhedralComplex; homogenized=false) =
  @view pm_object(PC).VERTICES[_vertex_indices(pm_object(PC)), (homogenized ? 1 : 2):end]

_matrix_for_polymake(::Val{_vertex_complex}) = _point_matrix

function _all_vertex_indices(P::Polymake.BigObject)
  vi = Polymake.get_attachment(P, "_all_vertex_indices")
  if isnothing(vi)
    A = P.VERTICES
    vi = Polymake.Vector{Polymake.to_cxx_type(Int64)}(Vector(1:Polymake.nrows(A)))
    Polymake.attach(P, "_all_vertex_indices", vi)
  end
  return vi
end

function _vertex_or_ray_complex(
  ::Type{Union{PointVector{T},RayVector{T}}}, PC::PolyhedralComplex{T}, i::Base.Integer
) where {T<:scalar_types}
  A = pm_object(PC).VERTICES
  if iszero(A[_all_vertex_indices(pm_object(PC))[i], 1])
    return ray_vector(
      coefficient_field(PC),
      @view pm_object(PC).VERTICES[_all_vertex_indices(pm_object(PC))[i], 2:end]
    )::RayVector{T}
  else
    return point_vector(
      coefficient_field(PC),
      @view pm_object(PC).VERTICES[_all_vertex_indices(pm_object(PC))[i], 2:end]
    )::PointVector{T}
  end
end

_vertices(
  as::Type{Union{RayVector{T},PointVector{T}}}, PC::PolyhedralComplex{T}
) where {T<:scalar_types} = SubObjectIterator{as}(
  PC, _vertex_or_ray_complex, length(_all_vertex_indices(pm_object(PC)))
)

@doc raw"""
    vertices_and_rays(PC::PolyhedralComplex)

Return the vertices and rays of `PC` as a combined set, up to lineality. This
function is mainly a helper function for [`maximal_polyhedra`](@ref maximal_polyhedra(PC::PolyhedralComplex{T}) where {T<:scalar_types}).

# Examples
```jldoctest
julia> IM = incidence_matrix([[1,2,3],[1,3,4]]);

julia> VR = [0 0 0; 1 0 0; 0 1 0; -1 0 0];

julia> far_vertices = [2,3,4];

julia> L = [0 0 1];

julia> PC = polyhedral_complex(IM, VR, far_vertices, L)
Polyhedral complex in ambient dimension 3

julia> for vr in vertices_and_rays(PC)
       println("$vr : $(typeof(vr))")
       end
QQFieldElem[0, 0, 0] : PointVector{QQFieldElem}
QQFieldElem[1, 0, 0] : RayVector{QQFieldElem}
QQFieldElem[0, 1, 0] : RayVector{QQFieldElem}
QQFieldElem[-1, 0, 0] : RayVector{QQFieldElem}
```
"""
vertices_and_rays(PC::PolyhedralComplex{T}) where {T<:scalar_types} =
  _vertices(Union{PointVector{T},RayVector{T}}, PC)

_vector_matrix(::Val{_vertex_or_ray_complex}, PC::PolyhedralComplex; homogenized=false) =
  homogenized ? pm_object(PC).VERTICES : @view pm_object(PC).VERTICES[:, 2:end]

_vertices(::Type{PointVector}, PC::PolyhedralComplex{T}) where {T<:scalar_types} =
  _vertices(PointVector{T}, PC)

vertices(PC::PolyhedralComplex{T}) where {T<:scalar_types} = vertices(PointVector{T}, PC)
_vertices(PC::PolyhedralComplex{T}) where {T<:scalar_types} = _vertices(PointVector{T}, PC)

_rays(as::Type{RayVector{T}}, PC::PolyhedralComplex{T}) where {T<:scalar_types} =
  SubObjectIterator{as}(
    PC, _ray_polyhedral_complex, length(_ray_indices_polyhedral_complex(pm_object(PC)))
  )

rays(as::Type{RayVector{T}}, PC::PolyhedralComplex{T}) where {T<:scalar_types} =
  lineality_dim(PC) == 0 ? _rays(RayVector{T}, PC) : _empty_subobjectiterator(as, PC)

@doc raw"""
    rays_modulo_lineality(as, PC::PolyhedralComplex)

Return the rays of the recession fan of `PC` up to lineality as a `NamedTuple`
with two iterators. If `PC` has lineality `L`, then the iterator
`rays_modulo_lineality` iterates over representatives of the rays of `PC/L`.
The iterator `lineality_basis` gives a basis of the lineality space `L`.

See also [`rays`](@ref rays(PC::PolyhedralComplex{T}) where {T<:scalar_types}) and [`lineality_space`](@ref lineality_space(PC::PolyhedralComplex{T}) where {T<:scalar_types}).

# Examples
```jldoctest
julia> VR = [0 0 0; 1 0 0; 0 1 0; -1 0 0];

julia> IM = incidence_matrix([[1,2,3],[1,3,4]]);

julia> far_vertices = [2,3,4];

julia> L = [0 0 1];

julia> PC = polyhedral_complex(IM, VR, far_vertices, L)
Polyhedral complex in ambient dimension 3

julia> RML = rays_modulo_lineality(PC)
(rays_modulo_lineality = RayVector{QQFieldElem}[[1, 0, 0], [0, 1, 0], [-1, 0, 0]], lineality_basis = RayVector{QQFieldElem}[[0, 0, 1]])

julia> RML.rays_modulo_lineality
3-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0, 0]
 [0, 1, 0]
 [-1, 0, 0]

julia> RML.lineality_basis
1-element SubObjectIterator{RayVector{QQFieldElem}}:
 [0, 0, 1]
```
"""
rays_modulo_lineality(PC::PolyhedralComplex{T}) where {T<:scalar_types} =
  rays_modulo_lineality(
    NamedTuple{
      (:rays_modulo_lineality, :lineality_basis),
      Tuple{SubObjectIterator{RayVector{T}},SubObjectIterator{RayVector{T}}},
    },
    PC,
  )

function rays_modulo_lineality(
  ::Type{
    NamedTuple{
      (:rays_modulo_lineality, :lineality_basis),
      Tuple{SubObjectIterator{RayVector{T}},SubObjectIterator{RayVector{T}}},
    },
  },
  PC::PolyhedralComplex{T},
) where {T<:scalar_types}
  return (
    rays_modulo_lineality=_rays(RayVector{T}, PC), lineality_basis=lineality_space(PC)
  )
end

rays_modulo_lineality(
  U::Type{RayVector{T}}, PC::PolyhedralComplex{T}
) where {T<:scalar_types} = _rays(U, PC)

@doc raw"""
    minimal_faces(as, PC::PolyhedralComplex)

Return the minimal faces of a polyhedral complex as a `NamedTuple` with two
iterators. For a polyhedral complex without lineality, the `base_points` are
the vertices. If `PC` has lineality `L`, then every minimal face is an affine
translation `p+L`, where `p` is only unique modulo `L`. The return type is a
dict, the key `:base_points` gives an iterator over such `p`, and the key
`:lineality_basis` lets one access a basis for the lineality space `L` of `PC`.

See also [`vertices`](@ref vertices(as::Type{PointVector{T}}, PC::PolyhedralComplex{T}) where {T<:scalar_types}) and [`lineality_space`](@ref).

# Examples
```jldoctest
julia> VR = [0 0 0; 1 0 0; 0 1 0; -1 0 0];

julia> IM = incidence_matrix([[1,2,3],[1,3,4]]);

julia> far_vertices = [2,3,4];

julia> L = [0 0 1];

julia> PC = polyhedral_complex(IM, VR, far_vertices, L)
Polyhedral complex in ambient dimension 3

julia> MFPC = minimal_faces(PC)
(base_points = PointVector{QQFieldElem}[[0, 0, 0]], lineality_basis = RayVector{QQFieldElem}[[0, 0, 1]])

julia> MFPC.base_points
1-element SubObjectIterator{PointVector{QQFieldElem}}:
 [0, 0, 0]

julia> MFPC.lineality_basis
1-element SubObjectIterator{RayVector{QQFieldElem}}:
 [0, 0, 1]
```
"""
minimal_faces(PC::PolyhedralComplex{T}) where {T<:scalar_types} = minimal_faces(
  NamedTuple{
    (:base_points, :lineality_basis),
    Tuple{SubObjectIterator{PointVector{T}},SubObjectIterator{RayVector{T}}},
  },
  PC,
)

function minimal_faces(
  ::Type{
    NamedTuple{
      (:base_points, :lineality_basis),
      Tuple{SubObjectIterator{PointVector{T}},SubObjectIterator{RayVector{T}}},
    },
  },
  PC::PolyhedralComplex{T},
) where {T<:scalar_types}
  return (base_points=_vertices(PointVector{T}, PC), lineality_basis=lineality_space(PC))
end

minimal_faces(as::Type{PointVector{T}}, PC::PolyhedralComplex{T}) where {T<:scalar_types} =
  _vertices(PointVector{T}, PC)

@doc raw"""
    rays([as::Type{T} = RayVector,] PC::PolyhedralComplex)

Return the rays of `PC`. The rays are defined to be the far vertices, i.e. the
one-dimensional faces of the recession cones of its polyhedra, so if `PC` has
lineality, there are no rays.

See also [`rays_modulo_lineality`](@ref rays_modulo_lineality(PC::PolyhedralComplex{T}) where {T<:scalar_types}) and [`vertices`](@ref vertices(as::Type{PointVector{T}}, PC::PolyhedralComplex{T}) where {T<:scalar_types}).

Optional arguments for `as` include
* `RayVector`.

# Examples
```jldoctest
julia> IM = incidence_matrix([[1,2,3],[1,3,4]]);

julia> VR = [0 0; 1 0; 1 1; 0 1];

julia> PC = polyhedral_complex(IM, VR, [2])
Polyhedral complex in ambient dimension 2

julia> rays(PC)
1-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0]

julia> matrix(QQ, rays(RayVector, PC))
[1   0]
```
The following complex has no vertices:
```
julia> IM = incidence_matrix([[1,2,3],[1,3,4]]);

julia> VR = [0 0 0; 1 0 0; 0 1 0; -1 0 0];

julia> far_vertices = [2,3,4];

julia> L = [0 0 1];

julia> PC = polyhedral_complex(IM, VR, far_vertices, L)
Polyhedral complex in ambient dimension 3

julia> rays(PC)
0-element SubObjectIterator{RayVector{QQFieldElem}}
```
"""
rays(PC::PolyhedralComplex{T}) where {T<:scalar_types} = rays(RayVector{T}, PC)
_rays(PC::PolyhedralComplex{T}) where {T<:scalar_types} = _rays(RayVector{T}, PC)

rays(as::Type{<:RayVector}, PC::PolyhedralComplex{T}) where {T<:scalar_types} =
  rays(RayVector{T}, PC)

_ray_indices_polyhedral_complex(PC::Polymake.BigObject) =
  collect(Polymake.to_one_based_indexing(PC.FAR_VERTICES))

_ray_polyhedral_complex(
  U::Type{RayVector{T}}, PC::PolyhedralComplex{T}, i::Base.Integer
) where {T<:scalar_types} = ray_vector(
  coefficient_field(PC),
  @view pm_object(PC).VERTICES[_ray_indices_polyhedral_complex(pm_object(PC))[i], 2:end]
)::U

_matrix_for_polymake(::Val{_ray_polyhedral_complex}) = _vector_matrix

_vector_matrix(::Val{_ray_polyhedral_complex}, PC::PolyhedralComplex; homogenized=false) =
  @view pm_object(PC).VERTICES[
    _ray_indices_polyhedral_complex(pm_object(PC)), (homogenized ? 1 : 2):end
  ]

_maximal_polyhedron(
  ::Type{Polyhedron{T}}, PC::PolyhedralComplex{T}, i::Base.Integer
) where {T<:scalar_types} =
  Polyhedron{T}(Polymake.fan.polytope(pm_object(PC), i - 1), coefficient_field(PC))

_vertex_indices(::Val{_maximal_polyhedron}, PC::PolyhedralComplex) =
  pm_object(PC).MAXIMAL_POLYTOPES[:, _vertex_indices(pm_object(PC))]

_ray_indices(::Val{_maximal_polyhedron}, PC::PolyhedralComplex) =
  pm_object(PC).MAXIMAL_POLYTOPES[:, _ray_indices_polyhedral_complex(pm_object(PC))]

_incidencematrix(::Val{_maximal_polyhedron}) = _vertex_and_ray_indices

_vertex_and_ray_indices(::Val{_maximal_polyhedron}, PC::PolyhedralComplex) =
  pm_object(PC).MAXIMAL_POLYTOPES

@doc raw"""
    maximal_polyhedra(PC::PolyhedralComplex)

Return the maximal polyhedra of `PC`

Optionally `IncidenceMatrix` can be passed as a first argument to return the
incidence matrix specifying the maximal polyhedra of `PC`. The indices returned
refer to the output of [`vertices_and_rays`](@ref vertices_and_rays(PC::PolyhedralComplex{T}) where {T<:scalar_types}).

# Examples
```jldoctest
julia> IM = incidence_matrix([[1,2,3],[1,3,4]])
2×4 IncidenceMatrix
[1, 2, 3]
[1, 3, 4]


julia> VR = [0 0; 1 0; 1 1; 0 1]
4×2 Matrix{Int64}:
 0  0
 1  0
 1  1
 0  1

julia> PC = polyhedral_complex(IM, VR, [2])
Polyhedral complex in ambient dimension 2

julia> maximal_polyhedra(PC)
2-element SubObjectIterator{Polyhedron{QQFieldElem}}:
 Polyhedron in ambient dimension 2
 Polytope in ambient dimension 2

julia> maximal_polyhedra(IncidenceMatrix, PC)
2×4 IncidenceMatrix
[1, 2, 3]
[1, 3, 4]
```
"""
maximal_polyhedra(PC::PolyhedralComplex{T}) where {T<:scalar_types} =
  SubObjectIterator{Polyhedron{T}}(PC, _maximal_polyhedron, n_maximal_polyhedra(PC))

@doc raw"""
    n_maximal_polyhedra(PC::PolyhedralComplex)

Return the number of maximal polyhedra of `PC`

# Examples
```jldoctest
julia> IM = incidence_matrix([[1,2,3],[1,3,4]])
2×4 IncidenceMatrix
[1, 2, 3]
[1, 3, 4]


julia> VR = [0 0; 1 0; 1 1; 0 1]
4×2 Matrix{Int64}:
 0  0
 1  0
 1  1
 0  1

julia> PC = polyhedral_complex(IM, VR, [2])
Polyhedral complex in ambient dimension 2

julia> n_maximal_polyhedra(PC)
2
```
"""
n_maximal_polyhedra(PC::PolyhedralComplex) = pm_object(PC).N_MAXIMAL_POLYTOPES

@doc raw"""
    is_simplicial(PC::PolyhedralComplex)

Determine whether the polyhedral complex is simplicial.

# Examples
```jldoctest
julia> IM = incidence_matrix([[1,2,3],[1,3,4]]);

julia> VR = [0 0; 1 0; 1 1; 0 1];

julia> PC = polyhedral_complex(IM, VR)
Polyhedral complex in ambient dimension 2

julia> is_simplicial(PC)
true
```
"""
is_simplicial(PC::PolyhedralComplex) = pm_object(PC).SIMPLICIAL::Bool

@doc raw"""
    is_pure(PC::PolyhedralComplex)

Determine whether the polyhedral complex is pure.

# Examples
```jldoctest
julia> IM = incidence_matrix([[1,2,3],[1,3,4]]);

julia> VR = [0 0; 1 0; 1 1; 0 1];

julia> PC = polyhedral_complex(IM, VR)
Polyhedral complex in ambient dimension 2

julia> is_pure(PC)
true
```
"""
is_pure(PC::PolyhedralComplex) = pm_object(PC).PURE::Bool

@doc raw"""
    dim(PC::PolyhedralComplex)

Compute the dimension of the polyhedral complex.

# Examples
```jldoctest
julia> IM = incidence_matrix([[1,2,3],[1,3,4]]);

julia> VR = [0 0; 1 0; 1 1; 0 1];

julia> PC = polyhedral_complex(IM, VR)
Polyhedral complex in ambient dimension 2

julia> dim(PC)
2
```
"""
dim(PC::PolyhedralComplex) = Polymake.fan.dim(pm_object(PC))::Int

@doc raw"""
    polyhedra_of_dim(PC::PolyhedralComplex, polyhedron_dim::Int)

Return the polyhedra of a given dimension in the polyhedral complex `PC`.

# Examples
```jldoctest
julia> IM = incidence_matrix([[1,2,3],[1,3,4]]);

julia> VR = [0 0; 1 0; 1 1; 0 1];

julia> PC = polyhedral_complex(IM, VR);

julia> P1s = polyhedra_of_dim(PC,1)
5-element SubObjectIterator{Polyhedron{QQFieldElem}}:
 Polytope in ambient dimension 2
 Polytope in ambient dimension 2
 Polytope in ambient dimension 2
 Polytope in ambient dimension 2
 Polytope in ambient dimension 2

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
function polyhedra_of_dim(
  PC::PolyhedralComplex{T}, polyhedron_dim::Int
) where {T<:scalar_types}
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
  return SubObjectIterator{Polyhedron{T}}(
    PC, _ith_polyhedron, length(rfaces), (f_dim=n, f_ind=rfaces)
  )
end

function _ith_polyhedron(
  ::Type{Polyhedron{T}},
  PC::PolyhedralComplex{T},
  i::Base.Integer;
  f_dim::Int=-1,
  f_ind::Vector{Int64}=Vector{Int64}(),
) where {T<:scalar_types}
  pface = Polymake.row(Polymake.fan.cones_of_dim(pm_object(PC), f_dim), f_ind[i])
  V = pm_object(PC).VERTICES[collect(pface), :]
  L = pm_object(PC).LINEALITY_SPACE
  PT = _scalar_type_to_polymake(T)
  return Polyhedron{T}(
    Polymake.polytope.Polytope{PT}(; VERTICES=V, LINEALITY_SPACE=L), coefficient_field(PC)
  )
end

@doc raw"""
    lineality_space(PC::PolyhedralComplex)

Return the lineality space of `PC`.

# Examples
```jldoctest
julia> VR = [0 0 0; 1 0 0; 0 1 0; -1 0 0];

julia> IM = incidence_matrix([[1,2,3],[1,3,4]]);

julia> far_vertices = [2,3,4];

julia> L = [0 0 1];

julia> PC = polyhedral_complex(IM, VR, far_vertices, L)
Polyhedral complex in ambient dimension 3

julia> lineality_space(PC)
1-element SubObjectIterator{RayVector{QQFieldElem}}:
 [0, 0, 1]
```
"""
lineality_space(PC::PolyhedralComplex{T}) where {T<:scalar_types} =
  SubObjectIterator{RayVector{T}}(PC, _lineality_complex, lineality_dim(PC))

_lineality_complex(
  U::Type{RayVector{T}}, PC::PolyhedralComplex{T}, i::Base.Integer
) where {T<:scalar_types} =
  ray_vector(coefficient_field(PC), @view pm_object(PC).LINEALITY_SPACE[i, 2:end])::U

_generator_matrix(::Val{_lineality_complex}, PC::PolyhedralComplex; homogenized=false) =
  if homogenized
    pm_object(PC).LINEALITY_SPACE
  else
    @view pm_object(PC).LINEALITY_SPACE[:, 2:end]
  end

_matrix_for_polymake(::Val{_lineality_complex}) = _generator_matrix

@doc raw"""
    lineality_dim(PC::PolyhedralComplex)

Return the lineality dimension of `PC`.

# Examples
```jldoctest
julia> VR = [0 0 0; 1 0 0; 0 1 0; -1 0 0];

julia> IM = incidence_matrix([[1,2,3],[1,3,4]]);

julia> far_vertices = [2,3,4];

julia> L = [0 0 1];

julia> PC = polyhedral_complex(IM, VR, far_vertices, L)
Polyhedral complex in ambient dimension 3

julia> lineality_dim(PC)
1
```
"""
lineality_dim(PC::PolyhedralComplex) = pm_object(PC).LINEALITY_DIM::Int

@doc raw"""
    f_vector(PC::PolyhedralComplex)

Compute the vector $(f₀,f₁,f₂,...,f_{dim(PC)})$ where $f_i$ is the number of
faces of $PC$ of dimension $i$.

# Examples
```jldoctest
julia> VR = [0 0; 1 0; -1 0; 0 1];

julia> IM = incidence_matrix([[1,2,4],[1,3,4]]);

julia> far_vertices = [2,3,4];

julia> PC = polyhedral_complex(IM, VR, far_vertices);

julia> f_vector(PC)
3-element Vector{Int64}:
 1
 3
 2
```
"""
function f_vector(PC::PolyhedralComplex)
  ldim = lineality_dim(PC)
  f_vec = vcat(zeros(Int64, ldim), [length(polyhedra_of_dim(PC, i)) for i in ldim:dim(PC)])
  return f_vec
end

@doc raw"""
    n_rays(PC::PolyhedralComplex)

Return the number of rays of `PC`.

# Examples
```jldoctest
julia> VR = [0 0; 1 0; -1 0; 0 1];

julia> IM = incidence_matrix([[1,2,4],[1,3,4]]);

julia> far_vertices = [2,3,4];

julia> PC = polyhedral_complex(IM, VR, far_vertices);

julia> n_rays(PC)
3
```
"""
n_rays(PC::PolyhedralComplex) = lineality_dim(PC) == 0 ? _n_rays(PC) : 0
_n_rays(PC::PolyhedralComplex) = length(pm_object(PC).FAR_VERTICES)

@doc raw"""
    n_vertices(PC::PolyhedralComplex)

Return the number of vertices of `PC`.

# Examples
```jldoctest
julia> VR = [0 0; 1 0; -1 0; 0 1];

julia> IM = incidence_matrix([[1,2,4],[1,3,4]]);

julia> far_vertices = [2,3,4];

julia> PC = polyhedral_complex(IM, VR, far_vertices);

julia> n_vertices(PC)
1
```
"""
n_vertices(PC::PolyhedralComplex) = lineality_dim(PC) == 0 ? _n_vertices(PC) : 0
_n_vertices(PC::PolyhedralComplex) = pm_object(PC).N_VERTICES - _n_rays(PC)

@doc raw"""
    n_polyhedra(PC::PolyhedralComplex)

Return the total number of polyhedra in the polyhedral complex `PC`.

# Examples
```jldoctest
julia> VR = [0 0; 1 0; -1 0; 0 1];

julia> IM = incidence_matrix([[1,2,4],[1,3,4]]);

julia> far_vertices = [2,3,4];

julia> PC = polyhedral_complex(IM, VR, far_vertices);

julia> n_polyhedra(PC)
6
```
"""
n_polyhedra(PC::PolyhedralComplex) = sum(f_vector(PC))

@doc raw"""
    codim(PC::PolyhedralComplex)

Compute the codimension of a polyhedral complex.

# Examples
```
julia> VR = [0 0; 1 0; -1 0; 0 1];

julia> IM = incidence_matrix([[1,2],[1,3],[1,4]]);

julia> far_vertices = [2,3,4];

julia> PC = polyhedral_complex(IM, VR, far_vertices)
A polyhedral complex in ambient dimension 2

julia> codim(PC)
1
```
"""
codim(PC::PolyhedralComplex) = ambient_dim(PC) - dim(PC)

@doc raw"""
    is_embedded(PC::PolyhedralComplex)

Return `true` if `PC` is embedded, i.e. if its vertices can be computed as a
subset of some $\mathbb{R}^n$.

# Examples
```jldoctest
julia> VR = [0 0; 1 0; -1 0; 0 1];

julia> IM = incidence_matrix([[1,2],[1,3],[1,4]]);

julia> PC = polyhedral_complex(IM, VR)
Polyhedral complex in ambient dimension 2

julia> is_embedded(PC)
true
```
"""
function is_embedded(PC::PolyhedralComplex)
  pmo = pm_object(PC)
  schedule = Polymake.call_method(pmo, :get_schedule, "VERTICES")
  return !isnothing(schedule)
end
