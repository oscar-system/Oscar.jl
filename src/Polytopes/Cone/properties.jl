###############################################################################
###############################################################################
### Iterators
###############################################################################
###############################################################################

struct ConeRayIterator
    cone::Cone
end

function Base.iterate(iter::ConeRayIterator, index = 1)
    rays = pm_cone(iter.cone).RAYS
    if size(rays, 1) < index
        return nothing
    end

    return (rays[index, :], index + 1)
end
Base.eltype(::Type{ConeRayIterator}) = Polymake.Vector{Polymake.Rational}
Base.length(iter::ConeRayIterator) = nrays(iter.cone)

"""
   rays(C::Cone)

Returns the rays of a cone.
"""
rays(C::Cone) = ConeRayIterator(C)

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

Returns the number of rays of the cone `C`.
"""
nrays(C::Cone) = pm_cone(C).N_RAYS

"""
   dim(C::Cone)

Returns the dimension of a cone.
"""
dim(C::Cone) = pm_cone(C).CONE_DIM

"""
   ambient_dim(C::Cone)

Returns the ambient dimension of a cone.
"""
ambient_dim(C::Cone) = pm_cone(C).CONE_AMBIENT_DIM

"""
   codim(C::Cone)

Returns the codimension of a cone.
"""
codim(C::Cone) = ambient_dim(C)-dim(C)

###############################################################################
## Boolean properties
###############################################################################
"""
   ispointed(C::Cone)

   Determine whether the cone is pointed, i.e. whether 0 is a face of the cone
"""
ispointed(C::Cone) = pm_cone(C).POINTED

"""
   isfulldimensional(C::Cone)

   Determine whether the cone is full dimensional
"""
isfulldimensional(C::Cone) = pm_cone(C).FULL_DIM

###############################################################################
## Points properties
###############################################################################

"""
   rays_as_point_matrix(C::Cone)

   Return the rays of a cone as rows in a matrix.
"""
function rays_as_point_matrix(C::Cone)
    pm_cone(C).RAYS
end


"""
   facets_as_point_matrix(C::Cone)

Returns the facets of a cone as rows of a matrix.
"""
facets_as_point_matrix(C::Cone) = pm_cone(C).FACETS


"""
   lineality_space(C::Cone)

   Returns a basis of the lineality space of a cone.
"""
lineality_space(C::Cone) = pm_cone(C).LINEALITY_SPACE

"""
   hilbert_basis(C::Cone)

   Returns the Hilbert basis of a pointed cone as the rows of a matrix.
"""
function hilbert_basis(C::Cone)
   if ispointed(C)
      return pm_cone(C).HILBERT_BASIS_GENERATORS[1]
   else
      throw(ArgumentError("Cone not pointed."))
   end
end
