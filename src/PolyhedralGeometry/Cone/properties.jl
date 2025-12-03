###############################################################################
###############################################################################
### Iterators
###############################################################################
###############################################################################

rays(as::Type{RayVector{T}}, C::Cone{T}) where {T<:scalar_types} =
  lineality_dim(C) == 0 ? _rays(as, C) : _empty_subobjectiterator(as, C)
_rays(as::Type{RayVector{T}}, C::Cone{T}) where {T<:scalar_types} =
  SubObjectIterator{as}(C, _ray_cone, _n_rays(C))

_ray_cone(U::Type{RayVector{T}}, C::Cone{T}, i::Base.Integer) where {T<:scalar_types} =
  ray_vector(coefficient_field(C), view(pm_object(C).RAYS, i, :))::U

_vector_matrix(::Val{_ray_cone}, C::Cone; homogenized=false) =
  homogenized ? homogenize(coefficient_field(C), pm_object(C).RAYS, 0) : pm_object(C).RAYS

_matrix_for_polymake(::Val{_ray_cone}) = _vector_matrix

rays(::Type{<:RayVector}, C::Cone{T}) where {T<:scalar_types} = rays(RayVector{T}, C)
_rays(::Type{<:RayVector}, C::Cone{T}) where {T<:scalar_types} = _rays(RayVector{T}, C)

@doc raw"""
    rays([as::Type{T} = RayVector,] C::Cone)

Return the rays of `C` in the format defined by `as`. The rays are defined to be the
one-dimensional faces, so if `C` has lineality, there are no rays.

See also [`rays_modulo_lineality`](@ref rays_modulo_lineality(C::Cone{T}) where {T<:scalar_types}).

Optional arguments for `as` include
* `RayVector`.

# Examples
Here a cone is constructed from three rays. Calling `rays` reveals that one of
these was redundant:
```jldoctest
julia> R = [1 0; 0 1; 0 2];

julia> PO = positive_hull(R);

julia> rays(PO)
2-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0]
 [0, 1]
```

The rays can also be converted to a matrix using the `matrix(ring, ...)` function.
If `ring=ZZ` the primitive generators of the rays are returned.
```jldoctest
julia> R = [1 0; 2 3];

julia> P = positive_hull(R);

julia> rays(P)
2-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0]
 [1, 3//2]

julia> matrix(QQ, rays(RayVector, P))
[1      0]
[1   3//2]

julia> matrix(ZZ, rays(P))
[1   0]
[2   3]
```
A half-space has no rays:
```
julia> UH = cone_from_inequalities([-1 0 0])
Polyhedral cone in ambient dimension 3

julia> rays(UH)
0-element SubObjectIterator{RayVector{QQFieldElem}}
```
"""
rays(C::Cone{T}) where {T<:scalar_types} = rays(RayVector{T}, C)
_rays(C::Cone{T}) where {T<:scalar_types} = _rays(RayVector{T}, C)

@doc raw"""                                                 
    rays_modulo_lineality(as, C::Cone)

Return the rays of the cone of `C` up to lineality as a `NamedTuple` with two
iterators. If `C` has lineality `L`, then the iterator `rays_modulo_lineality`
iterates over representatives of the rays of `C/L`. The iterator
`lineality_basis` gives a basis of the lineality space `L`.

See also [`rays`](@ref rays(C::Cone{T}) where {T<:scalar_types}) and [`lineality_space`](@ref lineality_space(C::Cone{T}) where {T<:scalar_types}).

# Examples
For a pointed cone, with two generators, we get the usual rays:
```jldoctest
julia> C = positive_hull([1 0; 0 1])
Polyhedral cone in ambient dimension 2

julia> rays(C)
2-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0]
 [0, 1]

julia> RML = rays_modulo_lineality(C)
(rays_modulo_lineality = RayVector{QQFieldElem}[[1, 0], [0, 1]], lineality_basis = RayVector{QQFieldElem}[])

julia> RML.rays_modulo_lineality
2-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0]
 [0, 1]

julia> RML.lineality_basis
0-element SubObjectIterator{RayVector{QQFieldElem}}
```
If the cone has lineality, the second iterator iterates over a basis for the
space of lineality.  The following example has one generator for the positive hull plus one generator for the lineality space: 
```jldoctest
julia> C = positive_hull([1 0],[0 1])
Polyhedral cone in ambient dimension 2

julia> lineality_dim(C)
1

julia> rays(C)
0-element SubObjectIterator{RayVector{QQFieldElem}}

julia> RML = rays_modulo_lineality(C)
(rays_modulo_lineality = RayVector{QQFieldElem}[[1, 0]], lineality_basis = RayVector{QQFieldElem}[[0, 1]])

julia> RML.lineality_basis
1-element SubObjectIterator{RayVector{QQFieldElem}}:
 [0, 1]
```
"""
rays_modulo_lineality(C::Cone{T}) where {T<:scalar_types} = rays_modulo_lineality(
  NamedTuple{
    (:rays_modulo_lineality, :lineality_basis),
    Tuple{SubObjectIterator{RayVector{T}},SubObjectIterator{RayVector{T}}},
  },
  C,
)
function rays_modulo_lineality(
  ::Type{
    NamedTuple{
      (:rays_modulo_lineality, :lineality_basis),
      Tuple{SubObjectIterator{RayVector{T}},SubObjectIterator{RayVector{T}}},
    },
  },
  C::Cone{T},
) where {T<:scalar_types}
  return (rays_modulo_lineality=_rays(C), lineality_basis=lineality_space(C))
end
rays_modulo_lineality(::Type{<:RayVector}, C::Cone) = _rays(C)

@doc raw"""
    faces(C::Cone, face_dim::Int)

Return an iterator over the faces of `C` of dimension `face_dim`.

# Examples
Each 2-dimensional face of the 3-dimensional positive orthant is generated by
two pairwise distinct unit vectors.
```jldoctest
julia> PO = cone_from_inequalities([-1 0 0; 0 -1 0; 0 0 -1])
Polyhedral cone in ambient dimension 3

julia> for f in faces(PO, 2)
       println(rays(f))
       end
RayVector{QQFieldElem}[[0, 1, 0], [0, 0, 1]]
RayVector{QQFieldElem}[[1, 0, 0], [0, 0, 1]]
RayVector{QQFieldElem}[[1, 0, 0], [0, 1, 0]]
```
"""
function faces(C::Cone{T}, face_dim::Int) where {T<:scalar_types}
  face_dim == dim(C) - 1 &&
    return SubObjectIterator{Cone{T}}(C, _face_cone_facet, n_facets(C))
  n = face_dim - length(lineality_space(C))
  n < 1 && return nothing
  return SubObjectIterator{Cone{T}}(
    C, _face_cone, size(Polymake.polytope.faces_of_dim(pm_object(C), n), 1), (f_dim=n,)
  )
end

function _face_cone(
  ::Type{Cone{T}}, C::Cone{T}, i::Base.Integer; f_dim::Int=0
) where {T<:scalar_types}
  R = pm_object(C).RAYS[
    collect(
      Polymake.to_one_based_indexing(Polymake.polytope.faces_of_dim(pm_object(C), f_dim)[i])
    ),
    :,
  ]
  L = pm_object(C).LINEALITY_SPACE
  PT = _scalar_type_to_polymake(T)
  return Cone{T}(
    Polymake.polytope.Cone{PT}(; RAYS=R, LINEALITY_SPACE=L), coefficient_field(C)
  )
end

function _ray_indices(::Val{_face_cone}, C::Cone; f_dim::Int=0)
  f = Polymake.to_one_based_indexing(Polymake.polytope.faces_of_dim(pm_object(C), f_dim))
  return IncidenceMatrix([collect(f[i]) for i in 1:length(f)])
end

function _face_cone_facet(
  ::Type{Cone{T}}, C::Cone{T}, i::Base.Integer
) where {T<:scalar_types}
  R = pm_object(C).RAYS[collect(pm_object(C).RAYS_IN_FACETS[i, :]), :]
  L = pm_object(C).LINEALITY_SPACE
  PT = _scalar_type_to_polymake(T)
  return Cone{T}(
    Polymake.polytope.Cone{PT}(; RAYS=R, LINEALITY_SPACE=pm_object(C).LINEALITY_SPACE),
    coefficient_field(C),
  )
end

_ray_indices(::Val{_face_cone_facet}, C::Cone) = pm_object(C).RAYS_IN_FACETS

_incidencematrix(::Val{_face_cone}) = _ray_indices

_incidencematrix(::Val{_face_cone_facet}) = _ray_indices

###############################################################################
###############################################################################
### Access properties
###############################################################################
###############################################################################

###############################################################################
## Scalar properties
###############################################################################
@doc raw"""
    n_facets(C::Cone)

Return the number of facets of a cone `C`.

# Examples
The cone over a square at height one has four facets.
```jldoctest
julia> C = positive_hull([1 0 0; 1 1 0; 1 1 1; 1 0 1])
Polyhedral cone in ambient dimension 3

julia> n_facets(C)
4
```
"""
n_facets(C::Cone) = size(pm_object(C).FACETS, 1)::Int

@doc raw"""
    n_rays(C::Cone)

Return the number of rays of `C`.

# Examples
Here a cone is constructed from three rays. Calling `number_of_rays` reveals that one of these was redundant:
```jldoctest
julia> R = [1 0; 0 1; 0 2];

julia> PO = positive_hull(R);

julia> n_rays(PO)
2
```
"""
n_rays(C::Cone) = lineality_dim(C) == 0 ? _n_rays(C) : 0
_n_rays(C::Cone) = size(pm_object(C).RAYS, 1)::Int

@doc raw"""
    dim(C::Cone)

Return the dimension of `C`.

# Examples
The cone `C` in this example is 2-dimensional within a 3-dimensional ambient space.
```jldoctest
julia> C = positive_hull([1 0 0; 1 1 0; 0 1 0]);

julia> dim(C)
2
```
"""
dim(C::Cone) = pm_object(C).CONE_DIM::Int

@doc raw"""
    ambient_dim(C::Cone)

Return the ambient dimension of `C`.

# Examples
The cone `C` in this example is 2-dimensional within a 3-dimensional ambient space.
```jldoctest
julia> C = positive_hull([1 0 0; 1 1 0; 0 1 0]);

julia> ambient_dim(C)
3
```
"""
ambient_dim(C::Cone) = pm_object(C).CONE_AMBIENT_DIM::Int

@doc raw"""
    codim(C::Cone)

Return the codimension of `C`.

# Examples
The cone `C` in this example is 2-dimensional within a 3-dimensional ambient space.
```jldoctest
julia> C = positive_hull([1 0 0; 1 1 0; 0 1 0]);

julia> codim(C)
1
```
"""
codim(C::Cone) = ambient_dim(C) - dim(C)

@doc raw"""
    f_vector(C::Cone)

Compute the vector $(f₁,f₂,...,f_{dim(C) - 1})$ where $f_i$ is the number of
faces of $C$ of dimension $i$.

# Examples
Take the cone over a square, then the f-vector of the cone is the same as of the square.
```jldoctest
julia> C = positive_hull([1 0 0; 1 1 0; 1 1 1; 1 0 1])
Polyhedral cone in ambient dimension 3

julia> f_vector(C)
2-element Vector{ZZRingElem}:
 4
 4

julia> square = cube(2)
Polytope in ambient dimension 2

julia> f_vector(square)
2-element Vector{ZZRingElem}:
 4
 4
```
"""
function f_vector(C::Cone)
  pmc = pm_object(C)
  ldim = pmc.LINEALITY_DIM
  return Vector{ZZRingElem}(vcat(fill(0, ldim), pmc.F_VECTOR))
end

@doc raw"""
    lineality_dim(C::Cone)

Compute the dimension of the lineality space of $C$, i.e. the largest linear
subspace contained in $C$.

# Examples
A cone is pointed if and only if the dimension of its lineality space is zero.
```jldoctest
julia> C = positive_hull([1 0 0; 1 1 0; 1 1 1; 1 0 1])
Polyhedral cone in ambient dimension 3

julia> is_pointed(C)
true

julia> lineality_dim(C)
0

julia> C1 = positive_hull([1 0],[0 1; 0 -1])
Polyhedral cone in ambient dimension 2

julia> is_pointed(C1)
false

julia> lineality_dim(C1)
1
```
"""
lineality_dim(C::Cone) = pm_object(C).LINEALITY_DIM::Int

@doc raw"""
    facet_degrees(C::Cone)

Facet degrees of the cone. The degree of a facet is the number of adjacent facets. 
In particular a general $2$-dimensional cone has two facets (rays) that meet at the origin. 

# Examples
Produce the facet degrees of a cone over a square and a cone over a square pyramid. 
```jldoctest
julia> c = positive_hull([1 1 0; 1 -1 0; 1 0 1; 1 0 -1])
Polyhedral cone in ambient dimension 3

julia> facet_degrees(c)
4-element Vector{Int64}:
 2
 2
 2
 2

julia> c = positive_hull([1 0 1 0; 1 0 -1 0; 1 0 0 1; 1 0 0 -1; 1 1 0 0])
Polyhedral cone in ambient dimension 4

julia> facet_degrees(c)
5-element Vector{Int64}:
 4
 3
 3
 3
 3
```
"""
facet_degrees(C::Cone) = Vector{Int}(Polymake.polytope.facet_degrees(pm_object(C)))

@doc raw"""
    ray_degrees(C::Cone)

Ray degrees of the cone. If the cone has lineality, the output is empty since
there are no rays that are also faces. 

# Examples
```jldoctest
julia> c = cone_from_inequalities([-1 0 0; 0 -1 0])
Polyhedral cone in ambient dimension 3

julia> ray_degrees(c)
Int64[]

julia> c = positive_hull([1 0 1 0; 1 0 -1 0; 1 0 0 1; 1 0 0 -1; 1 1 0 0])
Polyhedral cone in ambient dimension 4

julia> ray_degrees(c) 
5-element Vector{Int64}:
 3
 3
 3
 3
 4
```
"""
function ray_degrees(C::Cone)
  pm_object(C).LINEALITY_DIM > 0 && return Vector{Int}()
  return Vector{Int}(Polymake.polytope.vertex_degrees(pm_object(C)))
end

###############################################################################
## Boolean properties
###############################################################################
@doc raw"""
    is_simplicial(C::Cone)

Determine whether `C` is simplicial, i.e. whether the number of ray generators
is the same as the dimension of the cone modulo lineality.

# Examples
```jldoctest
julia> C0 = positive_hull([0 1])
Polyhedral cone in ambient dimension 2

julia> is_simplicial(C0)
true

julia> C1 = positive_hull([1 0 0; 1 1 0; 1 1 1; 1 0 1])
Polyhedral cone in ambient dimension 3

julia> is_simplicial(C1)
false
```
"""
is_simplicial(C::Cone) = pm_object(C).SIMPLICIAL_CONE::Bool

@doc raw"""
    is_smooth(C::Cone{QQFieldElem})

Determine whether `C` is smooth, i.e. whether its ray generators form part of a
lattice basis.

# Examples
```jldoctest
julia> C0 = positive_hull([0 1])
Polyhedral cone in ambient dimension 2

julia> is_smooth(C0)
true

julia> C1 = positive_hull([1 1; 1 -1])
Polyhedral cone in ambient dimension 2

julia> is_smooth(C1)
false
```
"""
is_smooth(C::Cone{QQFieldElem}) = pm_object(C).SMOOTH_CONE::Bool

@doc raw"""
    is_pointed(C::Cone)

Determine whether `C` is pointed, i.e. whether the origin is a face of `C`.

# Examples
A cone with lineality is not pointed, but a cone only consisting of a single ray is.
```jldoctest
julia> C = positive_hull([1 0], [0 1]);

julia> is_pointed(C)
false

julia> C = positive_hull([1 0]);

julia> is_pointed(C)
true
```
"""
is_pointed(C::Cone) = pm_object(C).POINTED::Bool

@doc raw"""
    is_fulldimensional(C::Cone)

Determine whether `C` is full-dimensional.

# Examples
The cone `C` in this example is 2-dimensional within a 3-dimensional ambient space.
```jldoctest
julia> C = positive_hull([1 0 0; 1 1 0; 0 1 0]);

julia> is_fulldimensional(C)
false
```
"""
is_fulldimensional(C::Cone) = pm_object(C).FULL_DIM::Bool

###############################################################################
## Points properties
###############################################################################

# TODO: facets as `Vector`? or `Matrix`?
@doc raw"""
    facets(as::Type{T} = LinearHalfspace, C::Cone)

Return the facets of `C` in the format defined by `as`.

The allowed values for `as` are
* `Halfspace` (or its subtype `LinearHalfspace`),
* `Hyperplane` (or its subtype `LinearHyperplane1),
* `Cone`.

# Examples
```jldoctest
julia> c = positive_hull([1 0 0; 0 1 0; 1 1 1])
Polyhedral cone in ambient dimension 3

julia> f = facets(Halfspace, c)
3-element SubObjectIterator{LinearHalfspace{QQFieldElem}} over the halfspaces of R^3 described by:
-x_3 <= 0
-x_1 + x_3 <= 0
-x_2 + x_3 <= 0
```
"""
facets(
  as::Type{<:Union{LinearHalfspace{T},LinearHyperplane{T},Cone{T}}}, C::Cone{T}
) where {T<:scalar_types} =
  SubObjectIterator{as}(C, _facet_cone, n_facets(C))

_facet_cone(
  U::Type{LinearHalfspace{T}}, C::Cone{T}, i::Base.Integer
) where {T<:scalar_types} =
  linear_halfspace(coefficient_field(C), -pm_object(C).FACETS[[i], :])::U
_facet_cone(
  U::Type{LinearHyperplane{T}}, C::Cone{T}, i::Base.Integer
) where {T<:scalar_types} =
  linear_hyperplane(coefficient_field(C), -pm_object(C).FACETS[[i], :])::U

_facet_cone(::Type{Cone{T}}, C::Cone{T}, i::Base.Integer) where {T<:scalar_types} =
  Cone{T}(Polymake.polytope.facet(pm_object(C), i - 1), coefficient_field(C))

_linear_inequality_matrix(::Val{_facet_cone}, C::Cone) = -pm_object(C).FACETS

_linear_matrix_for_polymake(::Val{_facet_cone}) = _linear_inequality_matrix

_ray_indices(::Val{_facet_cone}, C::Cone) = pm_object(C).RAYS_IN_FACETS

_incidencematrix(::Val{_facet_cone}) = _ray_indices

facets(C::Cone{T}) where {T<:scalar_types} = facets(LinearHalfspace{T}, C)

facets(::Type{<:Halfspace}, C::Cone{T}) where {T<:scalar_types} =
  facets(LinearHalfspace{T}, C)
facets(::Type{<:Hyperplane}, C::Cone{T}) where {T<:scalar_types} =
  facets(LinearHyperplane{T}, C)

facets(::Type{Cone}, C::Cone{T}) where {T<:scalar_types} = facets(Cone{T}, C)

@doc raw"""
    lineality_space(C::Cone)

Return a basis of the lineality space of `C`.

# Examples
Three rays are used here to construct the upper half-plane. Actually, two of these rays point in opposite directions.
This gives us a 1-dimensional lineality.
```jldoctest
julia> UH = positive_hull([1 0; 0 1; -1 0]);

julia> lineality_space(UH)
1-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0]
```
"""
lineality_space(C::Cone{T}) where {T<:scalar_types} =
  SubObjectIterator{RayVector{T}}(C, _lineality_cone, lineality_dim(C))

_lineality_cone(
  U::Type{RayVector{T}}, C::Cone{T}, i::Base.Integer
) where {T<:scalar_types} =
  ray_vector(coefficient_field(C), view(pm_object(C).LINEALITY_SPACE, i, :))::U

_generator_matrix(::Val{_lineality_cone}, C::Cone; homogenized=false) =
  if homogenized
    homogenize(coefficient_field(C), pm_object(C).LINEALITY_SPACE, 0)
  else
    pm_object(C).LINEALITY_SPACE
  end

_matrix_for_polymake(::Val{_lineality_cone}) = _generator_matrix

@doc raw"""
    linear_span(C::Cone)

Return the (linear) hyperplanes generating the linear span of `C`.

# Examples
This 2-dimensional cone in $\mathbb{R}^3$ lives in exactly one hyperplane $H$, with
$H = \{ (x_1, x_2, x_3) | x_3 = 0 \}$.
```jldoctest
julia> c = positive_hull([1 0 0; 0 1 0]);

julia> linear_span(c)
1-element SubObjectIterator{LinearHyperplane{QQFieldElem}} over the hyperplanes of R^3 described by:
x_3 = 0
```
"""
linear_span(C::Cone{T}) where {T<:scalar_types} =
  SubObjectIterator{LinearHyperplane{T}}(C, _linear_span, size(pm_object(C).LINEAR_SPAN, 1))

_linear_span(
  U::Type{LinearHyperplane{T}}, C::Cone{T}, i::Base.Integer
) where {T<:scalar_types} =
  linear_hyperplane(coefficient_field(C), view(pm_object(C).LINEAR_SPAN, i, :))::U

_linear_equation_matrix(::Val{_linear_span}, C::Cone) = pm_object(C).LINEAR_SPAN

_linear_matrix_for_polymake(::Val{_linear_span}) = _linear_equation_matrix

@doc raw"""
    hilbert_basis(C::Cone{QQFieldElem})

Return the Hilbert basis of a pointed cone `C` as the rows of a matrix.

# Examples
This (non-smooth) cone in the plane has a hilbert basis with three elements.
```jldoctest; filter = r".*"
julia> C = positive_hull([1 0; 1 2])
A polyhedral cone in ambient dimension 2

julia> matrix(ZZ, hilbert_basis(C))
[1   0]
[1   2]
[1   1]

```
"""
function hilbert_basis(C::Cone{QQFieldElem})
  @req is_pointed(C) "Cone not pointed"
  return SubObjectIterator{PointVector{ZZRingElem}}(
    C, _hilbert_generator, size(pm_object(C).HILBERT_BASIS_GENERATORS[1], 1)
  )
end

_hilbert_generator(
  T::Type{PointVector{ZZRingElem}}, C::Cone{QQFieldElem}, i::Base.Integer
) = point_vector(ZZ, view(pm_object(C).HILBERT_BASIS_GENERATORS[1], i, :))::T

_generator_matrix(::Val{_hilbert_generator}, C::Cone; homogenized=false) =
  if homogenized
    homogenize(pm_object(C).HILBERT_BASIS_GENERATORS[1], 0)
  else
    pm_object(C).HILBERT_BASIS_GENERATORS[1]
  end

_matrix_for_polymake(::Val{_hilbert_generator}) = _generator_matrix

@doc raw"""
    issubset(C0::Cone, C1::Cone)                           
    
Check whether `C0` is a subset of the cone `C1`.
                                             
# Examples                 
```jldoctest                                                                                
julia> C0 = positive_hull([1 1])
Polyhedral cone in ambient dimension 2

julia> C1 = positive_hull([1 0; 0 1])
Polyhedral cone in ambient dimension 2

julia> issubset(C0, C1)
true

julia> issubset(C1, C0)
false
```
"""
Base.issubset(C0::Cone{T}, C1::Cone{T}) where {T<:scalar_types} =
  Polymake.polytope.included_polyhedra(pm_object(C0), pm_object(C1))::Bool

@doc raw"""
    in(v::AbstractVector, C::Cone)

Check whether the vector `v` is contained in the cone `C`.

# Examples
The positive orthant only contains vectors with non-negative entries:
```jldoctest
julia> C = positive_hull([1 0; 0 1]);

julia> [1, 2] in C
true

julia> [1, -2] in C
false
```
"""
Base.in(v::AbstractVector, C::Cone) = Polymake.polytope.contains(pm_object(C), v)::Bool

@doc raw"""
    relative_interior_point(C::Cone)

Compute a point in the relative interior point of `C`, i.e. a point in `C` not
contained in any facet.
"""
relative_interior_point(C::Cone{T}) where {T<:scalar_types} = point_vector(
  coefficient_field(C), view(Polymake.common.dense(pm_object(C).REL_INT_POINT), :)
)::PointVector{T} # broadcast_view
