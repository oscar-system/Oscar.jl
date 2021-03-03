###############################################################################
###############################################################################
### Definition and constructors
###############################################################################
###############################################################################

@doc Markdown.doc"""
    Cone(Rays)

     A polyhedral cone, not necessarily pointed, defined by the positive hull
     of the rays, Rays.
"""
struct Cone #a real polymake polyhedron
    pm_cone::Polymake.BigObjectAllocated
end
function Cone(Rays::Union{Oscar.MatElem,AbstractMatrix})
    Cone(Polymake.polytope.Cone{Polymake.Rational}(
        INPUT_RAYS = matrix_for_polymake(Rays),
    ))
end
function Cone(Rays::Union{Oscar.MatElem,AbstractMatrix}, LS::Union{Oscar.MatElem,AbstractMatrix}; non_redundant::Bool=false)
   if non_redundant
       Cone(Polymake.polytope.Cone{Polymake.Rational}(
           RAYS = matrix_for_polymake(Rays),
           LINEALITY_SPACE = matrix_for_polymake(LS),
       ))
   else
       Cone(Polymake.polytope.Cone{Polymake.Rational}(
           INPUT_RAYS = matrix_for_polymake(Rays),
           INPUT_LINEALITY = matrix_for_polymake(LS),
       ))
   end
end

==(C0::Cone, C1::Cone) = Polymake.polytope.equal_polyhedra(pm_cone(C0), pm_cone(C1))


"""
    positive_hull(generators::Union{Oscar.MatElem,AbstractMatrix})

A polyhedral cone, not necessarily pointed, defined by the positive hull
of the `generators`. Redundant rays are allowed in the generators.
"""
function positive_hull(generators::Union{Oscar.MatElem,AbstractMatrix})
    # TODO: Filter out zero rows
    C=Polymake.polytope.Cone{Polymake.Rational}(INPUT_RAYS =
      matrix_for_polymake(remove_zero_rows(generators)))
    Cone(C)
end


"""
    pm_cone(C::Cone)

Get the underlying polymake `Cone`.
"""
pm_cone(C::Cone) = C.pm_cone


###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################

function Base.show(io::IO, C::Cone)
    print(io,"A polyhedral cone of dimension $(dim(C))")
end

Polymake.visual(C::Cone; opts...) = Polymake.visual(pm_cone(C); opts...)

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
Base.eltype(::Type{ConeRayIterator}) = Polymake.VectorAllocated{Polymake.Rational}
Base.length(iter::ConeRayIterator) = n_rays(iter.cone)

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
   n_rays(C::Cone)

Returns the number of rays of the cone `C`.
"""
n_rays(C::Cone) = pm_cone(C).N_RAYS

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

