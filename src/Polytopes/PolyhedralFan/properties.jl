###############################################################################
###############################################################################
### Iterators
###############################################################################
###############################################################################

struct PolyhedralFanRayIterator
    fan::PolyhedralFan
end

function Base.iterate(iter::PolyhedralFanRayIterator, index = 1)
    rays = pm_fan(iter.fan).RAYS
    if size(rays, 1) < index
        return nothing
    end

    return (rays[index, :], index + 1)
end
Base.eltype(::Type{PolyhedralFanRayIterator}) = Polymake.Vector{Polymake.Rational}
Base.length(iter::PolyhedralFanRayIterator) = nrays(iter.fan)

"""
    rays(PF::PolyhedralFan)

Returns the rays of a polyhedral fan.
"""
rays(PF::PolyhedralFan) = PolyhedralFanRayIterator(PF)


"""
    maximal_cones(PF::PolyhedralFan)

Returns the maximal cones of a polyhedral fan.
"""

#TODO: should the documentation mention maximal_cones_as_incidence_matrix?
#      similarly for cone ray iterators and facet iterators?
@doc Markdown.doc"""
    maximal_cones(PF::PolyhedralFan, as = :cones)

Returns an iterator over the maximal cones of the polyhedral fan `PF`.
"""
function maximal_cones(PF::PolyhedralFan)
   MaximalConeIterator(PF)
end

struct MaximalConeIterator
    PF::PolyhedralFan
end

function Base.iterate(iter::MaximalConeIterator, index = 1)
    n_max_cones = nmaximal_cones(iter.PF)
    if index > n_max_cones
        return nothing
    end
    current_cone = Cone(Polymake.fan.cone(pm_fan(iter.PF), index - 1))
    return (current_cone, index + 1)
end
Base.length(iter::MaximalConeIterator) = nmaximal_cones(iter.PF)

###############################################################################
###############################################################################
### Access properties
###############################################################################
###############################################################################

###############################################################################
## Scalar properties
###############################################################################

"""
    dim(PF::PolyhedralFan)

Returns the dimension of a polyhedral fan.
"""
dim(PF::PolyhedralFan) = pm_fan(PF).FAN_DIM

"""
    nmaximal_cones(PF::PolyhedralFan)

Returns the number of maximal cones in a polyhedral fan `PF`.
"""
nmaximal_cones(PF::PolyhedralFan) = pm_fan(PF).N_MAXIMAL_CONES

"""
    ambient_dim(PF::PolyhedralFan)

Returns the ambient dimension of a polyhedral fan, which is the dimension of
the embedding space. This is equal to the dimension of the fan if the fan is
full-dimensional.
"""
ambient_dim(PF::PolyhedralFan) = pm_fan(PF).FAN_AMBIENT_DIM

"""
    nrays(PF::PolyhedralFan)

Returns the number of rays of a polyhedral fan.
"""
nrays(PF::PolyhedralFan) = pm_fan(PF).N_RAYS


###############################################################################
## Points properties
###############################################################################

"""
    lineality_space(PF::PolyhedralFan)

Returns the lineality_space of a polyhedral fan.
"""
lineality_space(PF::PolyhedralFan) = pm_fan(PF).LINEALITY_SPACE


"""
    rays_as_point_matrix(PF::PolyhedralFan)

Returns the rays of a polyhedral fan as rows of a matrix.
"""
rays_as_point_matrix(PF::PolyhedralFan) = pm_fan(PF).RAYS


"""
    maximal_cones_as_incidence_matrix(PF::PolyhedralFan)

Returns the maximal cones of a polyhedral fan as an incidence matrix where the
rows correspond to the maximal cones and the columns to the rays.
"""
function maximal_cones_as_incidence_matrix(PF::PolyhedralFan)
   IncidenceMatrix(pm_fan(PF).MAXIMAL_CONES)
end

###############################################################################
## Boolean properties
###############################################################################
"""
    issmooth(PF::PolyhedralFan)

Determine whether the fan is smooth.
"""
issmooth(PF::PolyhedralFan) = pm_fan(PF).SMOOTH_FAN

"""
    isregular(PF::PolyhedralFan)

Determine whether the fan is regular, i.e. the normal fan of a polytope.
"""
isregular(PF::PolyhedralFan) = pm_fan(PF).REGULAR

"""
    iscomplete(PF::PolyhedralFan)

Determine whether the fan is complete.
"""
iscomplete(PF::PolyhedralFan) = pm_fan(PF).COMPLETE
