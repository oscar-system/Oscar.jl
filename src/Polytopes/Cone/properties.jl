###############################################################################
###############################################################################
### Iterators
###############################################################################
###############################################################################

"""
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
rays(C::Cone) = VectorIterator{RayVector{Polymake.Rational}}(pm_cone(C).RAYS)

function faces(as::Type{T}, C::Cone, face_dim::Int) where T<:Union{Cone, Cones}
    if face_dim < 0
        return nothing
    end
    cfaces = Polymake.polytope.faces_of_dim(C.pm_cone,face_dim-length(lineality_space(C)))
    return PolyhedronOrConeIterator{AsTypeIdentities(as)}(C.pm_cone.RAYS,cfaces, C.pm_cone.LINEALITY_SPACE)
end

###############################################################################
###############################################################################
### Access properties
###############################################################################
###############################################################################

###############################################################################
## Scalar properties
###############################################################################

"""
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

"""
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

"""
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

"""
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

###############################################################################
## Boolean properties
###############################################################################
"""
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

"""
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

facets(C::Cone) = HalfspaceIterator{Halfspace}(pm_cone(C).FACETS)

"""
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

"""
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
      return VectorIterator{PointVector{Polymake.to_cxx_type(Int64)}}(pm_cone(C).HILBERT_BASIS_GENERATORS[1])
   else
      throw(ArgumentError("Cone not pointed."))
   end
end
