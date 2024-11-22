###############################################################################
###############################################################################
### Iterators
###############################################################################
###############################################################################

@doc raw"""
    rays([as::Type{T} = RayVector,] PF::PolyhedralFan)

Return the rays of `PF`. The rays are defined to be the
one-dimensional faces of its cones, so if `PF` has lineality, there are no rays.

See also [`rays_modulo_lineality`](@ref rays_modulo_lineality(F::_FanLikeType)).

Optional arguments for `as` include
* `RayVector`.

# Examples
The rays of a normal fan of a cube point in every positive and negative unit
direction.
```jldoctest
julia> C = cube(3);

julia> NF = normal_fan(C);

julia> rays(NF)
6-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0, 0]
 [-1, 0, 0]
 [0, 1, 0]
 [0, -1, 0]
 [0, 0, 1]
 [0, 0, -1]
```

As for the `Cone`, the rays may be converted to a matrix using the
`matrix(ring, ...)` function.
```jldoctest
julia> C = cube(3);

julia> NF = normal_fan(C);

julia> matrix(QQ, rays(NF))
[ 1    0    0]
[-1    0    0]
[ 0    1    0]
[ 0   -1    0]
[ 0    0    1]
[ 0    0   -1]
```
The following fan has no rays:
```
julia> IM = incidence_matrix([[1,2],[2,3]]);

julia> R = [1 0 0; 0 1 0; -1 0 0];

julia> L = [0 0 1];

julia> PF = polyhedral_fan(IM, R, L)
Polyhedral fan in ambient dimension 3

julia> rays(PF)
0-element SubObjectIterator{RayVector{QQFieldElem}}
```
"""
rays(PF::_FanLikeType) =
  if lineality_dim(PF) == 0
    _rays(PF)
  else
    _empty_subobjectiterator(RayVector{_get_scalar_type(PF)}, PF)
  end
_rays(PF::_FanLikeType) =
  SubObjectIterator{RayVector{_get_scalar_type(PF)}}(PF, _ray_fan, _n_rays(PF))

_ray_fan(U::Type{RayVector{T}}, PF::_FanLikeType, i::Base.Integer) where {T<:scalar_types} =
  ray_vector(coefficient_field(PF), view(pm_object(PF).RAYS, i, :))::U

_vector_matrix(::Val{_ray_fan}, PF::_FanLikeType; homogenized=false) =
  homogenized ? homogenize(pm_object(PF).RAYS, 0) : pm_object(PF).RAYS

_matrix_for_polymake(::Val{_ray_fan}) = _vector_matrix

_maximal_cone(::Type{Cone{T}}, PF::_FanLikeType, i::Base.Integer) where {T<:scalar_types} =
  Cone{T}(Polymake.fan.cone(pm_object(PF), i - 1), coefficient_field(PF))

@doc raw"""
    rays_modulo_lineality(as, F::PolyhedralFan)

Return the rays of the polyhedral fan `F` up to lineality as a `NamedTuple`
with two iterators. If `F` has lineality `L`, then the iterator
`rays_modulo_lineality` iterates over representatives of the rays of `F/L`.
The iterator `lineality_basis` gives a basis of the lineality space `L`.

See also [`rays`](@ref rays(PF::_FanLikeType)) and [`lineality_space`](@ref lineality_space(PF::_FanLikeType)).

# Examples
```jldoctest
julia> P = convex_hull(QQFieldElem, [0 0; 1 0])
Polyhedron in ambient dimension 2

julia> NF = normal_fan(P)
Polyhedral fan in ambient dimension 2

julia> rmlF = rays_modulo_lineality(NF)
(rays_modulo_lineality = RayVector{QQFieldElem}[[1, 0], [-1, 0]], lineality_basis = RayVector{QQFieldElem}[[0, 1]])

julia> rmlF.rays_modulo_lineality
2-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0]
 [-1, 0]

julia> rmlF.lineality_basis
1-element SubObjectIterator{RayVector{QQFieldElem}}:
 [0, 1]

julia> rays(NF)
0-element SubObjectIterator{RayVector{QQFieldElem}}
```
"""
rays_modulo_lineality(F::_FanLikeType) = rays_modulo_lineality(
  NamedTuple{
    (:rays_modulo_lineality, :lineality_basis),
    Tuple{
      SubObjectIterator{RayVector{_get_scalar_type(F)}},
      SubObjectIterator{RayVector{_get_scalar_type(F)}},
    },
  },
  F,
)

function rays_modulo_lineality(
  ::Type{
    NamedTuple{
      (:rays_modulo_lineality, :lineality_basis),
      Tuple{SubObjectIterator{RayVector{T}},SubObjectIterator{RayVector{T}}},
    },
  },
  F::_FanLikeType,
) where {T<:scalar_types}
  return (rays_modulo_lineality=_rays(F), lineality_basis=lineality_space(F))
end

rays_modulo_lineality(::Type{<:RayVector}, F::_FanLikeType) = _rays(F)

@doc raw"""
    maximal_cones(PF::PolyhedralFan)

Return the maximal cones of `PF`.

Optionally `IncidenceMatrix` can be passed as a first argument to return the
incidence matrix specifying the maximal cones of `PF`. In that case, the
indices refer to the output of [`rays_modulo_lineality(Cone)`](@ref rays_modulo_lineality(F::_FanLikeType)).

# Examples
Here we ask for the the number of rays for each maximal cone of the face fan of
the 3-cube and use that `maximal_cones` returns an iterator.
```jldoctest
julia> PF = face_fan(cube(3));

julia> for c in maximal_cones(PF)
         println(n_rays(c))
       end
4
4
4
4
4
4

julia> maximal_cones(IncidenceMatrix, PF)
6×8 IncidenceMatrix
[1, 3, 5, 7]
[2, 4, 6, 8]
[1, 2, 5, 6]
[3, 4, 7, 8]
[1, 2, 3, 4]
[5, 6, 7, 8]
```
"""
maximal_cones(PF::_FanLikeType) =
  SubObjectIterator{Cone{_get_scalar_type(PF)}}(PF, _maximal_cone, n_maximal_cones(PF))

_ray_indices(::Val{_maximal_cone}, PF::_FanLikeType) = pm_object(PF).MAXIMAL_CONES

_incidencematrix(::Val{_maximal_cone}) = _ray_indices

@doc raw"""
    cones(PF::PolyhedralFan, cone_dim::Int)

Return an iterator over the cones of `PF` of dimension `cone_dim`.

# Examples
The 12 edges of the 3-cube correspond to the 2-dimensional cones of its face fan:
```jldoctest
julia> PF = face_fan(cube(3));

julia> cones(PF, 2)
12-element SubObjectIterator{Cone{QQFieldElem}}:
 Polyhedral cone in ambient dimension 3
 Polyhedral cone in ambient dimension 3
 Polyhedral cone in ambient dimension 3
 Polyhedral cone in ambient dimension 3
 Polyhedral cone in ambient dimension 3
 Polyhedral cone in ambient dimension 3
 Polyhedral cone in ambient dimension 3
 Polyhedral cone in ambient dimension 3
 Polyhedral cone in ambient dimension 3
 Polyhedral cone in ambient dimension 3
 Polyhedral cone in ambient dimension 3
 Polyhedral cone in ambient dimension 3
```
"""
function cones(PF::_FanLikeType, cone_dim::Int)
  l = cone_dim - length(lineality_space(PF))
  t = Cone{_get_scalar_type(PF)}
  (l < 0 || dim(PF) == -1) && return _empty_subobjectiterator(t, PF)

  if l == 0
    return SubObjectIterator{t}(
      PF,
      (_, _, _) -> positive_hull(
        coefficient_field(PF), zeros(Int, ambient_dim(PF)), lineality_space(PF)
      ),
      1,
      NamedTuple(),
    )
  end

  return SubObjectIterator{t}(
    PF, _cone_of_dim, size(Polymake.fan.cones_of_dim(pm_object(PF), l), 1), (c_dim=l,)
  )
end

function _cone_of_dim(
  ::Type{Cone{T}}, PF::_FanLikeType, i::Base.Integer; c_dim::Int=0
) where {T<:scalar_types}
  R = pm_object(PF).RAYS[
    collect(Polymake.row(Polymake.fan.cones_of_dim(pm_object(PF), c_dim), i)), :,
  ]
  L = pm_object(PF).LINEALITY_SPACE
  PT = _scalar_type_to_polymake(T)
  return Cone{T}(
    Polymake.polytope.Cone{PT}(; RAYS=R, LINEALITY_SPACE=L), coefficient_field(PF)
  )
end

_ray_indices(::Val{_cone_of_dim}, PF::_FanLikeType; c_dim::Int=0) =
  Polymake.fan.cones_of_dim(pm_object(PF), c_dim)

_incidencematrix(::Val{_cone_of_dim}) = _ray_indices

@doc raw"""
    cones(PF::PolyhedralFan)

Return the ray indices of all non-zero-dimensional
cones in a polyhedral fan.

# Examples
```jldoctest
julia> PF = face_fan(cube(2))
Polyhedral fan in ambient dimension 2

julia> cones(PF)
8×4 IncidenceMatrix
[1, 3]
[2, 4]
[1, 2]
[3, 4]
[1]
[3]
[2]
[4]
```
"""
function cones(PF::_FanLikeType)
  pmo = pm_object(PF)
  ncones = pmo.HASSE_DIAGRAM.N_NODES
  cones = [Polymake._get_entry(pmo.HASSE_DIAGRAM.FACES, i) for i in 0:(ncones - 1)]
  cones = filter(x -> !(-1 in x) && length(x) > 0, cones)
  cones = [Polymake.to_one_based_indexing(x) for x in cones]
  return IncidenceMatrix([Vector{Int}(x) for x in cones])
end

###############################################################################
###############################################################################
### Access properties
###############################################################################
###############################################################################

###############################################################################
## Scalar properties
###############################################################################

@doc raw"""
    dim(PF::PolyhedralFan)

Return the dimension of `PF`.

# Examples
This fan in the plane contains a 2-dimensional cone and is thus 2-dimensional
itself.
```jldoctest
julia> PF = polyhedral_fan(incidence_matrix([[1, 2], [3]]), [1 0; 0 1; -1 -1]);

julia> dim(PF)
2
```
"""
dim(PF::_FanLikeType) = pm_object(PF).FAN_DIM::Int

@doc raw"""
    n_maximal_cones(PF::PolyhedralFan)

Return the number of maximal cones of `PF`.

# Examples
The cones given in this construction are non-redundant. Thus there are two
maximal cones.
```jldoctest
julia> PF = polyhedral_fan(incidence_matrix([[1, 2], [3]]), [1 0; 0 1; -1 -1]);

julia> n_maximal_cones(PF)
2
```
"""
n_maximal_cones(PF::_FanLikeType) = pm_object(PF).N_MAXIMAL_CONES::Int

@doc raw"""
    n_cones(PF::PolyhedralFan)

Return the number of cones of `PF`.

# Examples
The cones given in this construction are non-redundant. There are six
cones in this fan.
```jldoctest
julia> PF = polyhedral_fan(incidence_matrix([[1, 2], [3]]), [1 0; 0 1; -1 -1])
Polyhedral fan in ambient dimension 2

julia> n_cones(PF)
4
```
"""
n_cones(PF::_FanLikeType) = nrows(cones(PF))

@doc raw"""
    ambient_dim(PF::PolyhedralFan)

Return the ambient dimension `PF`, which is the dimension of the embedding
space.

This is equal to the dimension of the fan if and only if the fan is
full-dimensional.

# Examples
The normal fan of the 4-cube is embedded in the same ambient space.
```jldoctest
julia> ambient_dim(normal_fan(cube(4)))
4
```
"""
ambient_dim(PF::_FanLikeType) = pm_object(PF).FAN_AMBIENT_DIM::Int

@doc raw"""
    n_rays(PF::PolyhedralFan)

Return the number of rays of `PF`.

# Examples
The 3-cube has 8 vertices. Accordingly, its face fan has 8 rays.
```jldoctest
julia> n_rays(face_fan(cube(3)))
8
```
"""
n_rays(PF::_FanLikeType) = lineality_dim(PF) == 0 ? _n_rays(PF) : 0
_n_rays(PF::_FanLikeType) = pm_object(PF).N_RAYS::Int

@doc raw"""
    f_vector(PF::PolyhedralFan)

Compute the vector $(f₁,f₂,...,f_{dim(PF)-1})$ where $f_i$ is the number of
faces of $PF$ of dimension $i$.

# Examples
The f-vector of the normal fan of a polytope is the reverse of the f-vector of
the polytope.
```jldoctest
julia> c = cube(3)
Polytope in ambient dimension 3

julia> f_vector(c)
3-element Vector{ZZRingElem}:
 8
 12
 6


julia> nfc = normal_fan(c)
Polyhedral fan in ambient dimension 3

julia> f_vector(nfc)
3-element Vector{ZZRingElem}:
 6
 12
 8
```
"""
function f_vector(PF::_FanLikeType)
  pmf = pm_object(PF)
  ldim = pmf.LINEALITY_DIM
  return Vector{ZZRingElem}(vcat(fill(0, ldim), pmf.F_VECTOR))
end

@doc raw"""
    lineality_dim(PF::PolyhedralFan)

Return the dimension of the lineality space of the polyhedral fan `PF`, i.e.
the dimension of the largest linear subspace.

# Examples
The dimension of the lineality space is zero if and only if the fan is pointed.
```jldoctest
julia> C = convex_hull([0 0; 1 0])
Polyhedron in ambient dimension 2

julia> is_fulldimensional(C)
false

julia> nf = normal_fan(C)
Polyhedral fan in ambient dimension 2

julia> is_pointed(nf)
false

julia> lineality_dim(nf)
1
```
"""
lineality_dim(PF::_FanLikeType) = pm_object(PF).LINEALITY_DIM::Int

###############################################################################
## Points properties
###############################################################################

@doc raw"""
    lineality_space(PF::PolyhedralFan)

Return a non-redundant matrix whose rows are generators of the lineality space
of `PF`.

# Examples
This fan consists of two cones, one containing all the points with $y ≤ 0$ and
one containing all the points with $y ≥ 0$. The fan's lineality is the common
lineality of these two cones, i.e. in $x$-direction.
```jldoctest
julia> PF = polyhedral_fan(incidence_matrix([[1, 2, 3], [3, 4, 1]]), [1 0; 0 1; -1 0; 0 -1])
Polyhedral fan in ambient dimension 2

julia> lineality_space(PF)
1-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0]
```
"""
lineality_space(PF::_FanLikeType) =
  SubObjectIterator{RayVector{_get_scalar_type(PF)}}(PF, _lineality_fan, lineality_dim(PF))

_lineality_fan(
  U::Type{RayVector{T}}, PF::_FanLikeType, i::Base.Integer
) where {T<:scalar_types} =
  ray_vector(coefficient_field(PF), view(pm_object(PF).LINEALITY_SPACE, i, :))::U

_generator_matrix(::Val{_lineality_fan}, PF::_FanLikeType; homogenized=false) =
  homogenized ? homogenize(pm_object(PF).LINEALITY_SPACE, 0) : pm_object(PF).LINEALITY_SPACE

_matrix_for_polymake(::Val{_lineality_fan}) = _generator_matrix

###############################################################################
## Boolean properties
###############################################################################
@doc raw"""
    is_pointed(PF::PolyhedralFan)

Determine whether `PF` is pointed, i.e. all its cones are pointed.

# Examples
The normal fan of a non-fulldimensional polytope is not pointed.
```jldoctest
julia> C = convex_hull([0 0; 1 0])
Polyhedron in ambient dimension 2

julia> is_fulldimensional(C)
false

julia> nf = normal_fan(C)
Polyhedral fan in ambient dimension 2

julia> is_pointed(nf)
false

julia> lineality_dim(nf)
1
```
"""
is_pointed(PF::_FanLikeType) = pm_object(PF).POINTED::Bool

@doc raw"""
    is_smooth(PF::PolyhedralFan{QQFieldElem})

Determine whether `PF` is smooth.

# Examples
Even though the cones of this fan cover the positive orthant together, one of
these und thus the whole fan is not smooth.
```jldoctest
julia> PF = polyhedral_fan(incidence_matrix([[1, 2], [2, 3]]), [0 1; 2 1; 1 0]);

julia> is_smooth(PF)
false
```
"""
is_smooth(PF::_FanLikeTypeQQ) = pm_object(PF).SMOOTH_FAN::Bool

@doc raw"""
    is_regular(PF::PolyhedralFan)

Determine whether `PF` is regular, i.e. the normal fan of a polytope.

# Examples
This fan is not complete and thus not regular.
```jldoctest
julia> PF = polyhedral_fan(incidence_matrix([[1, 2], [3]]), [1 0; 0 1; -1 -1]);

julia> is_regular(PF)
false
```
"""
is_regular(PF::_FanLikeType) = pm_object(PF).REGULAR::Bool

@doc raw"""
    is_pure(PF::PolyhedralFan)

Determine whether `PF` is pure, i.e. all maximal cones have the same dimension.

# Examples
```jldoctest
julia> PF = polyhedral_fan(incidence_matrix([[1, 2], [3]]), [1 0; 0 1; -1 -1]);

julia> is_pure(PF)
false
```
"""
is_pure(PF::_FanLikeType) = pm_object(PF).PURE::Bool

@doc raw"""
    is_fulldimensional(PF::PolyhedralFan)

Determine whether `PF` is fulldimensional, i.e. at least one maximal cone has maximal
dimension.

# Examples
```jldoctest
julia> PF = polyhedral_fan(incidence_matrix([[1, 2], [3]]), [1 0; 0 1; -1 -1]);

julia> is_fulldimensional(PF)
true
```
"""
is_fulldimensional(PF::_FanLikeType) = pm_object(PF).FULL_DIM::Bool

@doc raw"""
    is_complete(PF::PolyhedralFan)

Determine whether `PF` is complete, i.e. its support, the set-theoretic union
of its cones, covers the whole space.

# Examples
Normal fans of polytopes are complete.
```jldoctest
julia> is_complete(normal_fan(cube(3)))
true
```
"""
is_complete(PF::_FanLikeType) = pm_object(PF).COMPLETE::Bool

@doc raw"""
    is_simplicial(PF::PolyhedralFan)

Determine whether `PF` is simplicial, i.e. every cone should be generated by a
basis of the ambient space.

# Examples
The `normal_fan` of the cube is simplicial, while the `face_fan` is not.
```jldoctest
julia> is_simplicial(normal_fan(cube(3)))
true

julia> is_simplicial(face_fan(cube(3)))
false
```
"""
is_simplicial(PF::_FanLikeType) = pm_object(PF).SIMPLICIAL::Bool

###############################################################################
## Primitive collections
###############################################################################

@doc raw"""
    primitive_collections(PF::PolyhedralFan)

Return the primitive collections of a polyhedral fan.

# Examples
```jldoctest
julia> primitive_collections(normal_fan(simplex(3)))
1-element Vector{Set{Int64}}:
 Set([4, 2, 3, 1])
```
"""
function primitive_collections(PF::_FanLikeType)
  @req is_simplicial(PF) "PolyhedralFan must be simplicial."
  I = ray_indices(maximal_cones(PF))
  K = simplicial_complex(I)
  return minimal_nonfaces(K)
end
