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


###############################################################################
###############################################################################
### Access properties
###############################################################################
###############################################################################

function rays_as_point_matrix(C::Cone)
    pm_cone(C).RAYS
end

"""
   n_rays(C::Cone)

Returns the number of rays of the cone `C`.
"""
n_rays(C::Cone) = size(pm_cone(C).RAYS, 1)

"""
   rays(C::Cone)

Returns the rays of a cone.
"""
rays(C::Cone) = ConeRayIterator(C)

"""
   dim(C::Cone)

Returns the dimension of a cone.
"""
dim(C::Cone) = C.pm_cone.CONE_DIM

"""
   ambient_dim(C::Cone)

Returns the ambient dimension of a cone.
"""
ambient_dim(C::Cone) = C.pm_cone.CONE_AMBIENT_DIM

"""
   codim(C::Cone)

Returns the codimension of a cone.
"""
codim(C::Cone) = ambient_dim(C)-dim(C)

"""
   facets(C::Cone)

Returns the facets of a cone.
"""
facets(C::Cone) = C.pm_cone.facets


"""
   lineality_space(C::Cone)

   Returns a basis of the lineality space of a cone.
"""
lineality_space(C::Cone) = C.pm_cone.LINEALITY_SPACE

