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

@doc raw"""
    primitive_generator(r::AbstractVector{T}) where T<:RationalUnion -> Vector{ZZRingElem}

Given a vector `r` which is not a zero vector, returns the primitive generator of `r`, meaning the integer point on the ray generated by `r` that is closest to the origin.

# Examples
```jldoctest
julia> PF = polyhedral_fan(IncidenceMatrix([[1, 2]]), [[2, -1], [0, 1]])
Polyhedral fan in ambient dimension 2

julia> r = rays(PF)[1]
2-element RayVector{QQFieldElem}:
 1
 -1//2

julia> primitive_generator(r)
2-element Vector{ZZRingElem}:
 2
 -1
```
"""
function primitive_generator(r::AbstractVector{T}) where {T<:RationalUnion}
  return Vector{ZZRingElem}(first(primitive_generator_with_scaling_factor(r)))
end

@doc raw"""
    primitive_generator_with_scaling_factor(r::AbstractVector{T}) where T<:RationalUnion -> Tuple{Vector{ZZRingElem}, QQFieldElem}

Given a vector `r` which is not a zero vector, returns the primitive generator `s` of `r` together with the rational number `q` such that `qr = s`.

# Examples
```jldoctest
julia> PF = polyhedral_fan(IncidenceMatrix([[1, 2]]), [[2, -1], [0, 1]])
Polyhedral fan in ambient dimension 2

julia> r = rays(PF)[1]
2-element RayVector{QQFieldElem}:
 1
 -1//2

julia> primitive_generator_with_scaling_factor(r)
(ZZRingElem[2, -1], 2)
```
"""
function primitive_generator_with_scaling_factor(
  r::AbstractVector{T}
) where {T<:RationalUnion}
  first_scaling_factor = ZZ(lcm(denominator.(r)))
  result = ZZ.(first_scaling_factor * r)
  g = gcd(result)
  @req g > 0 "The vector `r` cannot be a zero vector"
  scaling_factor = QQ(first_scaling_factor, g)
  result = map(x -> div(x, g), result)
  return Tuple{Vector{ZZRingElem},QQFieldElem}((result, scaling_factor))
end

@doc raw"""
    minimal_supercone_indices(PF::PolyhedralFan, v::AbstractVector{<:RationalUnion})
    -> Set{Int64}

Given an point $v$ inside the support of the polyhedral fan $PF$, return
the ray indices of the unique cone $σ$ in $PF$ such that $v$ is in the
relative interior of $σ$.
Note that if $v$ is a zero vector, then $v$ is in the relative interior
of the zero cone.

The cone $σ$ can be constructed by

```
RPF = matrix(coefficient_field(PF), rays(PF))
isempty(result) && (sigma = positive_hull([], lineality_space(PF)))
!isempty(result) && (sigma = positive_hull(RPF[[result...], :], lineality_space(PF)))
```

where `result` is the output of this function.

# Examples
```jldoctest
julia> PF = normal_fan(Oscar.simplex(3))
Polyhedral fan in ambient dimension 3

julia> v = [1, 1, 0]
3-element Vector{Int64}:
 1
 1
 0

julia> minimal_supercone_indices(PF, v)
Set{Int64} with 2 elements:
  2
  1
```
"""
function minimal_supercone_indices(
  PF::PolyhedralFan, v::AbstractVector{<:RationalUnion}
)
  m_cone_indices = Oscar._get_maximal_cones_containing_vector(PF, v)
  @req !is_empty(m_cone_indices) "Point `v` should be in the support of the fan `PF`."
  IPF = maximal_cones(IncidenceMatrix, PF)
  m_cone_ray_indices = intersect([row(IPF, i) for i in m_cone_indices]...)
  isempty(m_cone_ray_indices) && return Set{Int64}()
  RPF = matrix(coefficient_field(PF), rays(PF))
  m_cone = positive_hull(RPF[[m_cone_ray_indices...], :], lineality_space(PF))
  PFm = facets(m_cone)
  zero_facets = findall(f -> f.a * v == [0], PFm)
  Rm = matrix(coefficient_field(m_cone), rays(m_cone))
  RIPF = IncidenceMatrix(m_cone.pm_cone.RAYS_IN_FACETS)
  some_list = [row(RIPF, i) for i in zero_facets]
  is_empty(some_list) && return m_cone_ray_indices
  ray_indices_in_m_cone = intersect(some_list...)
  function corresponding_ray_in_PF(i)
    r_i = rays_modulo_lineality(m_cone).rays_modulo_lineality[i]
    findfirst(isequal(r_i), rays_modulo_lineality(PF).rays_modulo_lineality)
  end
  ray_indices_in_PF = Set{Int64}(
    map(corresponding_ray_in_PF, collect(ray_indices_in_m_cone))
  )
  return ray_indices_in_PF
end

@doc raw"""
    minimal_supercone_coordinates(PF::PolyhedralFan, v::AbstractVector{<:RationalUnion})
    -> Vector{QQFieldElem}

Let $PF$ be a pointed polyhedral fan and let $u_1, \ldots, u_n$ be the
minimal generators of the rays of $PF$.
Let $\sigma$ be the unique cone in $PF$ such that $v$ is in the relative
interior of $v$.
Note that if $v$ is a zero vector, then $v$ is in the relative interior
of the zero cone.
This function returns a vector $(p_1, \ldots, p_n)$ of nonnegative
rational numbers such that both of the following hold:
  * the vector $v$ is equal to $p_1 u_1 + \ldots + p_n u_n$, and
  * if $u_i$ is not in $\sigma$, then $p_i = 0$.

If $PF$ is simplicial, then $(p_1, \ldots, p_n)$ is unique.

# Examples
```jldoctest
julia> PF = normal_fan(Oscar.simplex(3))
Polyhedral fan in ambient dimension 3

julia> v = [1, 1, 0]
3-element Vector{Int64}:
 1
 1
 0

julia> minimal_supercone_indices(PF, v)
Set{Int64} with 2 elements:
  2
  1
```
"""
function minimal_supercone_coordinates(
  PF::PolyhedralFan, v::AbstractVector{<:RationalUnion}
)
  # This function probably only makes sense for fans with no lineality
  @req is_pointed(PF) "The polyhedral fan must be pointed."

  inds = sort(collect(minimal_supercone_indices(PF, v)))
  M_ZZ = matrix(ZZ, rays(PF))[inds, :]
  result = zeros(QQFieldElem, n_rays(PF))
  if is_simplicial(PF)
    M_QQ = matrix(QQ, M_ZZ)
    v_QQ = Vector{QQFieldElem}(QQ.(v))
    coords = solve(M_QQ, v_QQ)
  else
    # If the fan is not simplicial, then solving over QQ might not give
    # vector with nonnegative entries.
    # Therefore, we instead solve over the integers.
    # To guarantee the existence of a solution, we first need to scale
    # the vector `v` by the exponent of the torsion part of the quotient
    # group of $\mathbb{Z}^n$ by the subgroup generated by the minimal
    # generators of the rays of `PF`.
    if isempty(elementary_divisors(M_ZZ))
      e = 1
    else
      e = last(setdiff(elementary_divisors(M_ZZ), [0]))
    end
    v_with_same_denominator = Hecke.FakeFmpqMat(Vector{QQFieldElem}(QQFieldElem.(v)))
    v_scale_factor = e * denominator(v_with_same_denominator)
    v_scaled = transpose(e * numerator(v_with_same_denominator))
    C = identity_matrix(ZZ, length(inds))
    coords_scaled_matrix = solve_mixed(transpose(M_ZZ), v_scaled, C)
    coords_scaled = coords_scaled_matrix[1, :]
    coords = QQFieldElem(1//v_scale_factor) .* (QQFieldElem.(coords_scaled))
  end
  for i in 1:length(inds)
    result[inds[i]] = coords[i]
  end
  return result
end

@doc raw"""
    is_minimal_supercone_coordinate_vector(PF::PolyhedralFan, v::AbstractVector{<:RationalUnion})
    -> Bool

Given a pointed polyhedral fan `PF` and a vector `v` of length equal to
the number of rays of `PF`, this function checks that both of the
following are true:
  * all of the entries of `v` are nonnegative,
  * there exists a cone $\sigma$ in `PF` such that for every $i$, if the
    $i$-th entry of $v$ is positive, then the $i$-th ray of `PF` belongs
    to $\sigma$.

# Examples
```jldoctest
julia> PF = normal_fan(Oscar.simplex(3))
Polyhedral fan in ambient dimension 3

julia> is_minimal_supercone_coordinate_vector(PF, [1, 1, 1, 0])
true

julia> is_minimal_supercone_coordinate_vector(PF, [1, 1, 1, 1])
false
```
"""
function is_minimal_supercone_coordinate_vector(
  PF::PolyhedralFan, v::AbstractVector{<:RationalUnion}
)
  @req is_pointed(PF) "The polyhedral fan must be pointed."
  @assert length(v) == n_rays(PF) "Length of v must match the number of rays."
  isnothing(findfirst(x -> x < 0, v)) || return false
  positive_indices = [i for i in 1:length(v) if v[i] > 0]
  maximal_cones_incidence_matrix = maximal_cones(IncidenceMatrix, PF)
  maximal_cones_indices = map(
    i -> row(maximal_cones_incidence_matrix, i),
    1:n_rows(maximal_cones_incidence_matrix),
  )
  for i in 1:n_rows(maximal_cones_incidence_matrix)
    issubset(positive_indices, maximal_cones_indices[i]) && return true
  end
  return false
end

@doc raw"""
    standard_coordinates(PF::PolyhedralFan, v::AbstractVector{<:RationalUnion})
    -> Vector{QQFieldElem}

If `is_minimal_supercone_coordinate_vector(PF, v)` is true, then return the dot product of `v` and the vector of primitive generators of the rays of `PF`.

# Examples
```jldoctest
julia> PF = normal_fan(Oscar.simplex(3))
Polyhedral fan in ambient dimension 3

julia> standard_coordinates(PF, [1, 1, 0, 1])
3-element Vector{QQFieldElem}:
 0
 0
 -1
```
"""
function standard_coordinates(PF::PolyhedralFan, coords::AbstractVector{<:RationalUnion})
  @assert is_minimal_supercone_coordinate_vector(
    PF, coords
  ) "Input vector must be a minimal supercone coordinate vector"
  primitive_ray_generators = Vector{Vector{QQFieldElem}}(map(primitive_generator, rays(PF)))
  return Vector{QQFieldElem}(sum(coords .* primitive_ray_generators))
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
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
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
