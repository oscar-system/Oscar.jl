###############################################################################
###############################################################################
### Iterators
###############################################################################
###############################################################################

#TODO: take into account lineality space

@doc Markdown.doc"""
    faces(P::Polyhedron, face_dim::Int)

Return an iterator over the faces of `P` of dimension `face_dim`.

# Examples
A `Vector` containing the six sides of the 3-dimensional cube can be obtained
via the following input:
```jldoctest
julia> F = faces(cube(3), 2)
6-element SubObjectIterator{Polyhedron}:
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
```
"""
function faces(P::Polyhedron, face_dim::Int)
    n = face_dim - length(lineality_space(P))
    n < 0 && return nothing
    pfaces = Polymake.to_one_based_indexing(Polymake.polytope.faces_of_dim(pm_object(P), n))
    nfaces = length(pfaces)
    rfaces = Vector{Int64}(undef, nfaces - binomial(nrays(P), n + 1))
    nfarf = 0
    farf = Polymake.to_one_based_indexing(pm_object(P).FAR_FACE)
    for index in 1:nfaces
        if pfaces[index] <= farf
            nfarf += 1
        else
            rfaces[index - nfarf] = index
        end
    end
    return SubObjectIterator{Polyhedron}(pm_object(P), _face_polyhedron, length(rfaces), (f_dim = n, f_ind = rfaces))
end

function _face_polyhedron(::Type{Polyhedron}, P::Polymake.BigObject, i::Base.Integer; f_dim::Int = -1, f_ind::Vector{Int64} = Vector{Int64}())
    pface = Polymake.to_one_based_indexing(Polymake.polytope.faces_of_dim(P, f_dim))[f_ind[i]]
    return Polyhedron(Polymake.polytope.Polytope(VERTICES = P.VERTICES[collect(pface),:], LINEALITY_SPACE = P.LINEALITY_SPACE))
end

function _vertex_incidences(::Val{_face_polyhedron}, P::Polymake.BigObject; f_dim = -1, f_ind::Vector{Int64} = Vector{Int64}())
    return IncidenceMatrix(collect.(Polymake.to_one_based_indexing(Polymake.polytope.faces_of_dim(P, f_dim)[f_ind])))[:, _vertex_indices(P)]
end

function _ray_incidences(::Val{_face_polyhedron}, P::Polymake.BigObject; f_dim = -1, f_ind::Vector{Int64} = Vector{Int64}())
    return IncidenceMatrix(collect.(Polymake.to_one_based_indexing(Polymake.polytope.faces_of_dim(P, f_dim)[f_ind])[:, _ray_indices(P)]))
end

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

@doc Markdown.doc"""
    vertices(as, P)

Return an iterator over the vertices of `P` in the format defined by `as`.

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
vertices(as::Type{PointVector{T}}, P::Polyhedron) where T = SubObjectIterator{as}(pm_object(P), _vertex_polyhedron, length(_vertex_indices(pm_object(P))))

_vertex_polyhedron(::Type{PointVector{T}}, P::Polymake.BigObject, i::Base.Integer) where T = PointVector{T}(P.VERTICES[_vertex_indices(P)[i], 2:end])

_point_matrix(::Val{_vertex_polyhedron}, P::Polymake.BigObject) = P.VERTICES[_vertex_indices(P), 2:end]

_matrix_for_polymake(::Val{_vertex_polyhedron}) = _point_matrix

vertices(::Type{PointVector}, P::Polyhedron) = vertices(PointVector{Polymake.Rational}, P)

@doc Markdown.doc"""
    vertices(P::Polyhedron)

Return an iterator over the vertices of a polyhedron `P` as points.

# Examples
The following code computes the vertices of the Minkowski sum of a triangle and
a square:
```jldoctest
julia> P = simplex(2) + cube(2);

julia> vertices(P)
5-element SubObjectIterator{PointVector{Polymake.Rational}}:
 [-1, -1]
 [2, -1]
 [2, 1]
 [-1, 2]
 [1, 2]
```
"""
vertices(P::Polyhedron) = vertices(PointVector, P)


@doc Markdown.doc"""
    nrays(P::Polyhedron)

Return the number of rays of `P`.

# Examples
Reflecting the input, the upper half-plane indeed has one ray.
```jldoctest
julia> UH = convex_hull([0 0],[0 1],[1 0]);

julia> nrays(UH)
1
```
"""
nrays(P::Polyhedron) = length(pm_object(P).FAR_FACE)

@doc Markdown.doc"""
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
nvertices(P::Polyhedron) = pm_object(P).N_VERTICES - nrays(P)


@doc Markdown.doc"""
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
2-element SubObjectIterator{RayVector{Polymake.Rational}}:
 [1, 0]
 [0, 1]
```
"""
rays(as::Type{RayVector{T}}, P::Polyhedron) where T = SubObjectIterator{as}(pm_object(P), _ray_polyhedron, length(_ray_indices(pm_object(P))))

_ray_polyhedron(::Type{RayVector{T}}, P::Polymake.BigObject, i::Base.Integer) where T = RayVector{T}(P.VERTICES[_ray_indices(P)[i], 2:end])

_vector_matrix(::Val{_ray_polyhedron}, P::Polymake.BigObject) = P.VERTICES[_ray_indices(P), 2:end]

_matrix_for_polymake(::Val{_ray_polyhedron}) = _vector_matrix

rays(::Type{RayVector}, P::Polyhedron) = rays(RayVector{Polymake.Rational}, P)

@doc Markdown.doc"""
    rays(P::Polyhedron)

Return minimal set of generators of the cone of unbounded directions of `P`
(i.e. its rays) as points.

# Examples
We can verify that the positive orthant of the plane is generated by the two
rays in positive unit direction:
```jldoctest
julia> PO = convex_hull([0 0], [1 0; 0 1]);

julia> rays(PO)
2-element SubObjectIterator{RayVector{Polymake.Rational}}:
 [1, 0]
 [0, 1]
```
"""
rays(P::Polyhedron) = rays(RayVector,P)

@doc Markdown.doc"""
    nfacets(P::Polyhedron)

Return the number of facets of `P`.

# Examples
The number of facets of the 5-dimensional cross polytope can be retrieved via
the following line:
```jldoctest
julia> nfacets(cross(5))
32
```
"""
nfacets(P::Polyhedron) = pm_object(P).N_FACETS

@doc Markdown.doc"""
    facets(as::Type{T} = Halfspace, P::Polyhedron)

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
6-element SubObjectIterator{Polyhedron}:
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3

julia> facets(Halfspace, C)
6-element SubObjectIterator{AffineHalfspace}:
 The Halfspace of R^3 described by
1: -x₁ ≦ 1

 The Halfspace of R^3 described by
1: x₁ ≦ 1

 The Halfspace of R^3 described by
1: -x₂ ≦ 1

 The Halfspace of R^3 described by
1: x₂ ≦ 1

 The Halfspace of R^3 described by
1: -x₃ ≦ 1

 The Halfspace of R^3 described by
1: x₃ ≦ 1
```
"""
facets(as::Type{T}, P::Polyhedron) where {R, S, T<:Union{AffineHalfspace, Pair{R, S}, Polyhedron}} = SubObjectIterator{as}(pm_object(P), _facet_polyhedron, pm_object(P).N_FACETS)

function _facet_polyhedron(::Type{T}, C::Polymake.BigObject, i::Base.Integer) where {R, S, T<:Union{Polyhedron, AffineHalfspace, Pair{R, S}}}
    h = decompose_hdata(C.FACETS[[i], :])
    return T(h[1], h[2][])
end

_affine_inequality_matrix(::Val{_facet_polyhedron}, C::Polymake.BigObject) = -C.FACETS

_affine_matrix_for_polymake(::Val{_facet_polyhedron}) = _affine_inequality_matrix

_halfspace_matrix_pair(::Val{_facet_polyhedron}, C::Polymake.BigObject) = decompose_hdata(C.FACETS)

facets(::Type{Pair}, P::Polyhedron) = facets(Pair{Polymake.Matrix{Polymake.Rational}, Polymake.Rational}, P)

@doc Markdown.doc"""
    facets(P::Polyhedron)

Return the facets of `P` as halfspaces.

# Examples
We can retrieve the six facets of the 3-dimensional cube this way:
```jldoctest
julia> C = cube(3);

julia> facets(C)
6-element SubObjectIterator{AffineHalfspace}:
 The Halfspace of R^3 described by
1: -x₁ ≦ 1

 The Halfspace of R^3 described by
1: x₁ ≦ 1

 The Halfspace of R^3 described by
1: -x₂ ≦ 1

 The Halfspace of R^3 described by
1: x₂ ≦ 1

 The Halfspace of R^3 described by
1: -x₃ ≦ 1

 The Halfspace of R^3 described by
1: x₃ ≦ 1
```
"""
facets(P::Polyhedron) = facets(AffineHalfspace, P)

facets(::Type{Halfspace}, P::Polyhedron) = facets(AffineHalfspace, P)

###############################################################################
###############################################################################
### Access properties
###############################################################################
###############################################################################

###############################################################################
## Scalar properties
###############################################################################
@doc Markdown.doc"""
    lineality_dim(P::Polyhedron)

Return the dimension of the lineality space, i.e. the dimension of the largest
affine subspace contained in `P`.

# Examples
Polyhedron with one lineality direction.
```jldoctest
julia> C = convex_hull([0 0], [1 0], [1 1])
A polyhedron in ambient dimension 2

julia> lineality_dim(C)
1
```
"""
lineality_dim(P::Polyhedron) = pm_object(P).LINEALITY_DIM


@doc Markdown.doc"""
    volume(P::Polyhedron)

Return the (Euclidean) volume of `P`.

# Examples
```jldoctest
julia> C = cube(2);

julia> volume(C)
4
```
"""
volume(P::Polyhedron) = (pm_object(P)).VOLUME


@doc Markdown.doc"""
    lattice_volume(P::Polyhedron)

Return the lattice volume of `P`.

# Examples
```jldoctest
julia> C = cube(2);

julia> lattice_volume(C)
8
```
"""
lattice_volume(P::Polyhedron) = (pm_object(P)).LATTICE_VOLUME


@doc Markdown.doc"""
    normalized_volume(P::Polyhedron)

Return the (normalized) volume of `P`.

# Examples
```jldoctest
julia> C = cube(2);

julia> normalized_volume(C)
8
```
"""
normalized_volume(P::Polyhedron) = factorial(dim(P))*(pm_object(P)).VOLUME

@doc Markdown.doc"""
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
dim(P::Polyhedron) = Polymake.polytope.dim(pm_object(P))


@doc Markdown.doc"""
    lattice_points(P::Polyhedron)

Return the integer points contained in the bounded polyhedron `P`.

# Examples
```jldoctest
julia> S = 2 * simplex(2);

julia> lattice_points(S)
6-element SubObjectIterator{PointVector{Polymake.Integer}}:
 [0, 0]
 [0, 1]
 [0, 2]
 [1, 0]
 [1, 1]
 [2, 0]
```
"""
function lattice_points(P::Polyhedron)
    pm_object(P).BOUNDED || throw(ArgumentError("Polyhedron not bounded"))
    return SubObjectIterator{PointVector{Polymake.Integer}}(pm_object(P), _lattice_point, size(pm_object(P).LATTICE_POINTS_GENERATORS[1], 1))
end

_lattice_point(::Type{PointVector{Polymake.Integer}}, P::Polymake.BigObject, i::Base.Integer) = PointVector{Polymake.Integer}(P.LATTICE_POINTS_GENERATORS[1][i, 2:end])

_point_matrix(::Val{_lattice_point}, P::Polymake.BigObject) = P.LATTICE_POINTS_GENERATORS[1][:, 2:end]

_matrix_for_polymake(::Val{_lattice_point}) = _point_matrix

@doc Markdown.doc"""
    interior_lattice_points(P::Polyhedron)

Return the integer points contained in the interior of the bounded polyhedron
`P`.

# Examples
```jldoctest
julia> c = cube(3)
A polyhedron in ambient dimension 3

julia> interior_lattice_points(c)
1-element SubObjectIterator{PointVector{Polymake.Integer}}:
 [0, 0, 0]
```
"""
function interior_lattice_points(P::Polyhedron)
    pm_object(P).BOUNDED || throw(ArgumentError("Polyhedron not bounded"))
    return SubObjectIterator{PointVector{Polymake.Integer}}(pm_object(P), _interior_lattice_point, size(pm_object(P).INTERIOR_LATTICE_POINTS, 1))
end

_interior_lattice_point(::Type{PointVector{Polymake.Integer}}, P::Polymake.BigObject, i::Base.Integer) = PointVector{Polymake.Integer}(P.INTERIOR_LATTICE_POINTS[i, 2:end])

_point_matrix(::Val{_interior_lattice_point}, P::Polymake.BigObject) = P.INTERIOR_LATTICE_POINTS[:, 2:end]

_matrix_for_polymake(::Val{_interior_lattice_point}) = _point_matrix

@doc Markdown.doc"""
    boundary_lattice_points(P::Polyhedron)

Return the integer points contained in the boundary of the bounded polyhedron
`P`.

# Examples
```jldoctest
julia> c = polarize(cube(3))
A polyhedron in ambient dimension 3

julia> boundary_lattice_points(c)
6-element SubObjectIterator{PointVector{Polymake.Integer}}:
 [-1, 0, 0]
 [0, -1, 0]
 [0, 0, -1]
 [0, 0, 1]
 [0, 1, 0]
 [1, 0, 0]
```
"""
function boundary_lattice_points(P::Polyhedron)
    pm_object(P).BOUNDED || throw(ArgumentError("Polyhedron not bounded"))
    return SubObjectIterator{PointVector{Polymake.Integer}}(pm_object(P), _boundary_lattice_point, size(pm_object(P).BOUNDARY_LATTICE_POINTS, 1))
end

_boundary_lattice_point(::Type{PointVector{Polymake.Integer}}, P::Polymake.BigObject, i::Base.Integer) = PointVector{Polymake.Integer}(P.BOUNDARY_LATTICE_POINTS[i, 2:end])

_point_matrix(::Val{_boundary_lattice_point}, P::Polymake.BigObject) = P.BOUNDARY_LATTICE_POINTS[:, 2:end]

_matrix_for_polymake(::Val{_boundary_lattice_point}) = _point_matrix

@doc Markdown.doc"""
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
ambient_dim(P::Polyhedron) = Polymake.polytope.ambient_dim(pm_object(P))

@doc Markdown.doc"""
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
@doc Markdown.doc"""
    lineality_space(P::Polyhedron)

Return a matrix whose row span is the lineality space of `P`.

# Examples
Despite not being reflected in this construction of the upper half-plane,
its lineality in $x$-direction is recognized:
```jldoctest
julia> UH = convex_hull([0 0],[0 1; 1 0; -1 0]);

julia> lineality_space(UH)
1-element SubObjectIterator{RayVector{Polymake.Rational}}:
 [1, 0]
```
"""
lineality_space(P::Polyhedron) = SubObjectIterator{RayVector{Polymake.Rational}}(pm_object(P), _lineality_polyhedron, lineality_dim(P))

_lineality_polyhedron(::Type{RayVector{Polymake.Rational}}, P::Polymake.BigObject, i::Base.Integer) = RayVector(P.LINEALITY_SPACE[i, 2:end])

_generator_matrix(::Val{_lineality_polyhedron}, P::Polymake.BigObject) = P.LINEALITY_SPACE[:, 2:end]

_matrix_for_polymake(::Val{_lineality_polyhedron}) = _generator_matrix

@doc Markdown.doc"""
    affine_hull(P::Polytope)

Return the (affine) hyperplanes generating the affine hull of `P`.

# Examples
This triangle in $\mathbb{R}^4$ is contained in the plane defined by
$P = \{ (x_1, x_2, x_3, x_4) | x_3 = 2 ∧ x_4 = 5 \}$.
```jldoctest
julia> t = convex_hull([0 0 2 5; 1 0 2 5; 0 1 2 5]);

julia> affine_hull(t)
2-element SubObjectIterator{AffineHyperplane}:
 The Hyperplane of R^4 described by
1: x₃ = 2

 The Hyperplane of R^4 described by
1: x₄ = 5

```
"""
affine_hull(P::Polyhedron) = SubObjectIterator{AffineHyperplane}(pm_object(P), _affine_hull, size(pm_object(P).AFFINE_HULL, 1))

function _affine_hull(::Type{AffineHyperplane}, P::Polymake.BigObject, i::Base.Integer)
    h = decompose_hdata(-P.AFFINE_HULL[[i], :])
    return AffineHyperplane(h[1], h[2][])
end

_affine_equation_matrix(::Val{_affine_hull}, P::Polymake.BigObject) = -P.AFFINE_HULL

@doc Markdown.doc"""
    recession_cone(P::Polyhedron)

Return the recession cone of `P`.

# Examples
```jldoctest
julia> P = Polyhedron([1 -2; -1 1; -1 0; 0 -1],[2,1,1,1]);

julia> vertices(P)
3-element SubObjectIterator{PointVector{Polymake.Rational}}:
 [0, -1]
 [-1, 0]
 [-1, -1]

julia> recession_cone(P)
A polyhedral cone in ambient dimension 2

julia> rays(recession_cone(P))
2-element SubObjectIterator{RayVector{Polymake.Rational}}:
 [1, 1/2]
 [1, 1]
```
"""
recession_cone(P::Polyhedron) = Cone(Polymake.polytope.recession_cone(pm_object(P)))

###############################################################################
## Boolean properties
###############################################################################
@doc Markdown.doc"""
    isfeasible(P::Polyhedron)

Check whether `P` is feasible, i.e. non-empty.

# Examples
```jldoctest
julia> P = Polyhedron([1 -1; -1 1; -1 0; 0 -1],[-1,-1,1,1]);

julia> isfeasible(P)
false
```
"""
isfeasible(P::Polyhedron) = pm_object(P).FEASIBLE

@doc Markdown.doc"""
    contains(P::Polyhedron, v::AbstractVector)

Check whether `P` contains `v`.

# Examples
The positive orthant only contains vectors with non-negative entries:
```jldoctest
julia> PO = Polyhedron([-1 0; 0 -1], [0, 0]);

julia> contains(PO, [1, 2])
true

julia> contains(PO, [1, -2])
false
```
"""
contains(P::Polyhedron, v::AbstractVector) = Polymake.polytope.contains(pm_object(P), [1; v])

@doc Markdown.doc"""
    issmooth(P::Polyhedron)

Check whether `P` is smooth.

# Examples
A cube is always smooth.
```jldoctest
julia> C = cube(8);

julia> issmooth(C)
true
```
"""
issmooth(P::Polyhedron) = pm_object(P).SMOOTH


@doc Markdown.doc"""
    isnormal(P::Polyhedron)

Check whether `P` is normal.

# Examples
The 3-cube is normal.
```jldoctest
julia> C = cube(3)
A polyhedron in ambient dimension 3

julia> isnormal(C)
true
```
But this pyramid is not:
```jldoctest
julia> P = convex_hull([0 0 0; 0 1 1; 1 1 0; 1 0 1]);

julia> isnormal(P)
false
```
"""
isnormal(P::Polyhedron) = pm_object(P).NORMAL


@doc Markdown.doc"""
    isbounded(P::Polyhedron)

Check whether `P` is bounded.

# Examples
```jldoctest
julia> P = Polyhedron([1 -3; -1 1; -1 0; 0 -1],[1,1,1,1]);

julia> isbounded(P)
false
```
"""
isbounded(P::Polyhedron) = pm_object(P).BOUNDED


@doc Markdown.doc"""
    isfulldimensional(P::Polyhedron)

Check whether `P` is full-dimensional.

# Examples
```jldoctest
julia> V = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];

julia> isfulldimensional(convex_hull(V))
false
```
"""
isfulldimensional(P::Polyhedron) = pm_object(P).FULL_DIM

@doc Markdown.doc"""
    f_vector(P::Polyhedron)

Compute the vector $(f₀,f₁,f₂,...,f_{(dim(P)-1))$` where $f_i$ is the number of
faces of $P$ of dimension $i$.

# Examples
Here we compute the f-vector of the 5-cube:
```jldoctest
julia> f_vector(cube(5))
5-element Vector{Int64}:
 32
 80
 80
 40
 10
```
"""
function f_vector(P::Polyhedron)
    ldim = lineality_dim(P)
    f_vec=vcat(zeros(Int64, ldim), [length(faces(P,i)) for i in ldim:dim(P)-1])
    return f_vec
end


@doc Markdown.doc"""
    relative_interior_point(P::Polyhedron)

Compute a point in the relative interior point of `P`, i.e. a point in `P` not
contained in any facet.

# Examples
The square $[-1,1]^3$ has the origin as a relative interior point.
```jldoctest
julia> square = cube(2)
A polyhedron in ambient dimension 2

julia> relative_interior_point(square)
2-element PointVector{Polymake.Rational}:
 0
 0

julia> vertices(square)
4-element VectorIterator{PointVector{Polymake.Rational}}:
 [-1, -1]
 [1, -1]
 [-1, 1]
 [1, 1]
```
"""
relative_interior_point(P::Polyhedron) = PointVector(dehomogenize(Polymake.common.dense(pm_object(P).REL_INT_POINT)))

@doc Markdown.doc"""
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
function support_function(P::Polyhedron; convention = :max)
    function h(ω::AbstractVector)
        lp=LinearProgram(P,ω; convention = convention)
        return solve_lp(lp)[1]
    end
    return h
end

@doc Markdown.doc"""
    print_constraints(A::AnyVecOrMat, b::AbstractVector; trivial::Bool = false)

Pretty print the constraints given by $P(A,b) = \{ x |  Ax ≤ b \}$.

Trivial inequalities are counted but omitted. They are included if `trivial` is
set to `true`.

# Examples
The 3-cube is given by $-1 ≦ x_i ≦ 0 ∀ i ∈ \{1, 2, 3\}$.
```jldoctest
julia> print_constraints([-1 0 4 5; 4 4 4 3; 1 0 0 0; 0 0 0 0; 0 0 0 0; 9 9 9 9], [0, 1, 2, 3, -4, 5])
1: -x₁ + 4*x₃ + 5*x₄ ≦ 0
2: 4*x₁ + 4*x₂ + 4*x₃ + 3*x₄ ≦ 1
3: x₁ ≦ 2
5: 0 ≦ -4
6: 9*x₁ + 9*x₂ + 9*x₃ + 9*x₄ ≦ 5

julia> print_constraints([-1 0 4 5; 4 4 4 3; 1 0 0 0; 0 0 0 0; 0 0 0 0; 9 9 9 9], [0, 1, 2, 3, -4, 5]; trivial = true)
1: -x₁ + 4*x₃ + 5*x₄ ≦ 0
2: 4*x₁ + 4*x₂ + 4*x₃ + 3*x₄ ≦ 1
3: x₁ ≦ 2
4: 0 ≦ 3
5: 0 ≦ -4
6: 9*x₁ + 9*x₂ + 9*x₃ + 9*x₄ ≦ 5
```
"""
function print_constraints(A::AnyVecOrMat, b::AbstractVector; trivial::Bool = false, io::IO = stdout, cmp::String = "≦")
    for i in 1:length(b)
        terms = Vector{String}(undef, size(A)[2])
        first = true
        for j in 1:size(A)[2]
            if iszero(A[i, j])
                terms[j] = ""
            else
                if isone(abs(A[i, j]))
                    terms[j] = first ? string(isone(A[i, j]) ? "x" : "-x" , ['₀'+ d for d in digits(j)]...) :
                        string(isone(A[i, j]) ? " + x" : " - x", ['₀'+ d for d in digits(j)]...)
                else
                    terms[j] = first ? string(A[i, j], "*x", ['₀'+ d for d in digits(j)]...) :
                        string(A[i, j] < 0 ? " - " : " + ", abs(A[i, j]), "*x", ['₀'+ d for d in digits(j)]...)
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
        println(io, string(i, ": ", terms..., " ", cmp, " ", b[i]))
    end
end

@doc Markdown.doc"""
    print_constraints(P::Polyhedron; trivial::Bool = false)

Pretty print the constraints given by $P(A,b) = \{ x |  Ax ≤ b \}$.

Trivial inequalities are counted but omitted. They are included if `trivial` is
set to `true`.

# Examples
The 3-cube is given by $-1 ≦ x_i ≦ 0 ∀ i ∈ \{1, 2, 3\}$.
```jldoctest
julia> print_constraints(cube(3))
1: -x₁ ≦ 1
2: x₁ ≦ 1
3: -x₂ ≦ 1
4: x₂ ≦ 1
5: -x₃ ≦ 1
6: x₃ ≦ 1
```
"""
print_constraints(P::Polyhedron; trivial::Bool = false, io::IO = stdout) = print_constraints(halfspace_matrix_pair(facets(P))...; trivial = trivial, io = io)

print_constraints(H::Halfspace; trivial::Bool = false, io::IO = stdout) = print_constraints(vcat(H.a'), [negbias(H)]; trivial = trivial, io = io)

print_constraints(H::Hyperplane; trivial::Bool = false, io::IO = stdout) = print_constraints(vcat(H.a'), [negbias(H)]; trivial = trivial, io = io, cmp = "=")

function Base.show(io::IO, H::Halfspace)
    n = length(H.a)
    if iszero(H.a) && H.b >= 0
        print(io, "The trivial Halfspace, R^$n")
    else
        print(io, "The Halfspace of R^$n described by\n")
        print_constraints(H; io=io)
    end
end

function Base.show(io::IO, H::Hyperplane)
    n = length(H.a)
    b = negbias(H)
    if b == 0 && iszero(H.a)
        print(io, "The trivial Hyperplane, R^$n")
    else
        print(io, "The Hyperplane of R^$n described by\n")
        print_constraints(H; io = io)
    end
end
