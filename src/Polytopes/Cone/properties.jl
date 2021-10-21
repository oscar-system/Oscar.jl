###############################################################################
###############################################################################
### Iterators
###############################################################################
###############################################################################

rays(as::Type{RayVector{T}}, C::Cone) where T = VectorIterator{as}(C.pm_cone.RAYS)

rays(::Type{RayVector}, C::Cone) = rays(RayVector{Polymake.Rational}, C)

@doc Markdown.doc"""
    rays(C::Cone)

Return the rays of `C`.

# Examples
Here a cone is constructed from three rays. Calling `rays` reveals that one of these was redundant:
```jldoctest
julia> R = [1 0; 0 1; 0 2];

julia> PO = positive_hull(R);

julia> rays(PO)
2-element VectorIterator{RayVector{Polymake.Rational}}:
 [1, 0]
 [0, 1]
```
"""
rays(C::Cone) = rays(RayVector{Polymake.Rational}, C)

function faces(as::Type{T}, C::Cone, face_dim::Int) where T<:Cone
   n = face_dim-length(lineality_space(C))
   if n < 1
      return nothing
   end
   cfaces = Polymake.to_one_based_indexing(Polymake.polytope.faces_of_dim(C.pm_cone, n))
   return PolyhedronOrConeIterator{as}(C.pm_cone.RAYS,cfaces, C.pm_cone.LINEALITY_SPACE)
end

faces(C::Cone, face_dim::Int) = faces(Cone, C, face_dim)

###############################################################################
###############################################################################
### Access properties
###############################################################################
###############################################################################

###############################################################################
## Scalar properties
###############################################################################
@doc Markdown.doc"""
    nfacets(C::Cone)

Return the number of facets of a cone `C`.

# Examples
The cone over a square at height one has four facets.
```jldoctest
julia> C = positive_hull([1 0 0; 1 1 0; 1 1 1; 1 0 1])
A polyhedral cone in ambient dimension 3

julia> nfacets(C)
4
```
"""
nfacets(C::Cone) = pm_cone(C).N_FACETS


@doc Markdown.doc"""
    nrays(C::Cone)

Return the number of rays of `C`.

# Examples
Here a cone is constructed from three rays. Calling `nrays` reveals that one of these was redundant:
```jldoctest
julia> R = [1 0; 0 1; 0 2];

julia> PO = positive_hull(R);

julia> nrays(PO)
2
```
"""
nrays(C::Cone) = pm_cone(C).N_RAYS

@doc Markdown.doc"""
    dim(C::Cone)

Return the dimension of `C`.

# Examples
The cone `C` in this example is 2-dimensional within a 3-dimensional ambient space.
```jldoctest
julia> C = Cone([1 0 0; 1 1 0; 0 1 0]);

julia> dim(C)
2
```
"""
dim(C::Cone) = pm_cone(C).CONE_DIM

@doc Markdown.doc"""
    ambient_dim(C::Cone)

Return the ambient dimension of `C`.

# Examples
The cone `C` in this example is 2-dimensional within a 3-dimensional ambient space.
```jldoctest
julia> C = Cone([1 0 0; 1 1 0; 0 1 0]);

julia> ambient_dim(C)
3
```
"""
ambient_dim(C::Cone) = pm_cone(C).CONE_AMBIENT_DIM

@doc Markdown.doc"""
    codim(C::Cone)

Return the codimension of `C`.

# Examples
The cone `C` in this example is 2-dimensional within a 3-dimensional ambient space.
```jldoctest
julia> C = Cone([1 0 0; 1 1 0; 0 1 0]);

julia> codim(C)
1
```
"""
codim(C::Cone) = ambient_dim(C)-dim(C)


@doc Markdown.doc"""
    f_vector(C::Cone)

Compute the vector $(f₁,f₂,...,f_{(dim(C)-1))$` where $f_i$ is the number of
faces of $C$ of dimension $i$.

# Examples
Take the cone over a square, then the f-vector of the cone is the same as of the square.
```jldoctest
julia> C = positive_hull([1 0 0; 1 1 0; 1 1 1; 1 0 1])
A polyhedral cone in ambient dimension 3

julia> f_vector(C)
2-element Vector{Polymake.Integer}:
 4
 4

julia> square = cube(2)
A polyhedron in ambient dimension 2

julia> f_vector(square)
2-element Vector{Int64}:
 4
 4
```
"""
function f_vector(C::Cone)
    pmc = pm_cone(C)
    ldim = pmc.LINEALITY_DIM
    return vcat(fill(0,ldim),pmc.F_VECTOR)
end


@doc Markdown.doc"""
    lineality_dim(C::Cone)

Compute the dimension of the lineality space of $C$, i.e. the largest linear
subspace contained in $C$.

# Examples
A cone is pointed if and only if the dimension of its lineality space is zero.
```jldoctest
julia> C = positive_hull([1 0 0; 1 1 0; 1 1 1; 1 0 1])
A polyhedral cone in ambient dimension 3

julia> ispointed(C)
true

julia> lineality_dim(C)
0

julia> C1 = Cone([1 0],[0 1; 0 -1])
A polyhedral cone in ambient dimension 2

julia> ispointed(C1)
false

julia> lineality_dim(C1)
1
```
"""
lineality_dim(C::Cone) = pm_cone(C).LINEALITY_DIM



###############################################################################
## Boolean properties
###############################################################################
@doc Markdown.doc"""
    ispointed(C::Cone)

Determine whether `C` is pointed, i.e. whether the origin is a face of `C`.

# Examples
A cone with lineality is not pointed, but a cone only consisting of a single ray is.
```jldoctest
julia> C = Cone([1 0], [0 1]);

julia> ispointed(C)
false

julia> C = Cone([1 0]);

julia> ispointed(C)
true
```
"""
ispointed(C::Cone) = pm_cone(C).POINTED

@doc Markdown.doc"""
    isfulldimensional(C::Cone)

Determine whether `C` is full-dimensional.

# Examples
The cone `C` in this example is 2-dimensional within a 3-dimensional ambient space.
```jldoctest
julia> C = Cone([1 0 0; 1 1 0; 0 1 0]);

julia> isfulldimensional(C)
false
```
"""
isfulldimensional(C::Cone) = pm_cone(C).FULL_DIM

###############################################################################
## Points properties
###############################################################################

# TODO: facets as `Vector`? or `Matrix`?
@doc Markdown.doc"""
    facets(as::Type{T} = Halfspace, C::Cone)

Return the facets of `C` in the format defined by `as`.

The allowed values for `as` are
* `Halfspace`,
* `Cone.

# Examples
```jldoctest
julia> c = positive_hull([1 0 0; 0 1 0; 1 1 1])
A polyhedral cone in ambient dimension 3

julia> f = facets(Halfspace, c)
3-element HalfspaceIterator{Halfspace}:
 The Halfspace of R^3 described by
1: -x₃ ≦ 0

 The Halfspace of R^3 described by
1: -x₁ + x₃ ≦ 0

 The Halfspace of R^3 described by
1: -x₂ + x₃ ≦ 0
```
"""
facets(as::Type{T}, C::Cone) where T<:Union{Halfspace, Cone} = HalfspaceIterator{as}(-pm_cone(C).FACETS)

@doc Markdown.doc"""
    facets(as::Type{T} = Halfspace, C::Cone)

Return the facets of `C` in the format defined by `as`.

The allowed values for `as` are
* `Halfspace`,
* `Cone.

# Examples
```jldoctest
julia> c = positive_hull([1 0 0; 0 1 0; 1 1 1])
A polyhedral cone in ambient dimension 3

julia> f = facets(c)
3-element HalfspaceIterator{Halfspace}:
 The Halfspace of R^3 described by
1: -x₃ ≦ 0

 The Halfspace of R^3 described by
1: -x₁ + x₃ ≦ 0

 The Halfspace of R^3 described by
1: -x₂ + x₃ ≦ 0
```
"""
facets(C::Cone) = facets(Halfspace, C)

@doc Markdown.doc"""
    lineality_space(C::Cone)

Return a basis of the lineality space of `C`.

# Examples
Three rays are used here to construct the upper half-plane. Actually, two of these rays point in opposite directions.
This gives us a 1-dimensional lineality.
```jldoctest
julia> UH = Cone([1 0; 0 1; -1 0]);

julia> lineality_space(UH)
1-element VectorIterator{RayVector{Polymake.Rational}}:
 [1, 0]
```
"""
lineality_space(C::Cone) = VectorIterator{RayVector{Polymake.Rational}}(pm_cone(C).LINEALITY_SPACE)

@doc Markdown.doc"""
    linear_span(C::Cone)

Return the (linear) hyperplanes generating the linear span of `C`.

# Examples
This 2-dimensional cone in $\mathbb{R}^3$ lives in exactly one hyperplane $H$, with
$H = \{ (x_1, x_2, x_3) | x_3 = 0 \}$.
```jldoctest
julia> c = Cone([1 0 0; 0 1 0]);

julia> linear_span(c)
1-element HalfspaceIterator{Hyperplane}:
 The Hyperplane of R^3 described by
1: x₃ = 0

```
"""
linear_span(C::Cone) = HalfspaceIterator{Hyperplane}(C.pm_cone.LINEAR_SPAN)

@doc Markdown.doc"""
    hilbert_basis(C::Cone)

Return the Hilbert basis of a pointed cone `C` as the rows of a matrix.

# Examples
This (non-smooth) cone in the plane has a hilbert basis with three elements.
```jldoctest; filter = r".*"
julia> C = Cone([1 0; 1 2])
A polyhedral cone in ambient dimension 2

julia> hilbert_basis(C)
pm::Matrix<pm::Integer>
1 0
1 1
1 2
```
"""
function hilbert_basis(C::Cone)
   if ispointed(C)
      return VectorIterator{PointVector{Polymake.Integer}}(pm_cone(C).HILBERT_BASIS_GENERATORS[1])
   else
      throw(ArgumentError("Cone not pointed."))
   end
end
