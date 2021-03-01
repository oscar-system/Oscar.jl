@doc Markdown.doc"""
   PolyhedralFan(Rays, Cones)

    A polyhedral fan formed from the Rays by taking the cones corresponding to
    the indices of the sets in Cones.
"""
struct PolyhedralFan
   pm_fan::Polymake.BigObjectAllocated
end


###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, PF::PolyhedralFan)
    print(io, "A polyhedral fan in dimension $(ambient_dim(PF))")
end

function PolyhedralFan(Cones::Array{Cone,1})
   BigObjectArray = Polymake.Array{Polymake.BigObject}(length(Cones))
   for i in 1:length(Cones)
      BigObjectArray[i] = pm_cone(Cones[i])
   end
   PolyhedralFan(Polymake.fan.check_fan_objects(BigObjectArray))
end

#=
function PolyhedralFan(Rays::Union{Oscar.MatElem,AbstractMatrix}, Cones::AbstractArray{AbstractArray{Int64}})
   PolyhedralFan(Polymake.fan.PolyhedralFan{Polymake.Rational}(
      INPUT_RAYS = matrix_for_polymake(Rays),
      INPUT_CONES = [[x-1 for x in cone] for cone in Cones],
   ))
end
=#
function PolyhedralFan(Rays::Union{Oscar.MatElem,AbstractMatrix}, LS::Union{Oscar.MatElem,AbstractMatrix}, Incidence::IncidenceMatrix)
   PolyhedralFan(Polymake.fan.PolyhedralFan{Polymake.Rational}(
      INPUT_RAYS = matrix_for_polymake(Rays),
      INPUT_LINEALITY = matrix_for_polymake(LS),
      INPUT_CONES = [collect(to_zero_based_indexing(row(Incidence.pm_incidencematrix, i))) for i in 1:size(Incidence.pm_incidencematrix, 1)],
   ))
end

function PolyhedralFan(Rays::Union{Oscar.MatElem,AbstractMatrix}, Incidence::IncidenceMatrix)
   PolyhedralFan(Polymake.fan.PolyhedralFan{Polymake.Rational}(
      INPUT_RAYS = matrix_for_polymake(Rays),
      INPUT_CONES = [collect(to_zero_based_indexing(row(Incidence.pm_incidencematrix, i))) for i in 1:size(Incidence.pm_incidencematrix, 1)],
   ))
end

#Same construction for when the user gives Array{Bool,2} as incidence matrix
function PolyhedralFan(Rays::Union{Oscar.MatElem,AbstractMatrix}, LS::Union{Oscar.MatElem,AbstractMatrix}, Incidence::Array{Bool,2})
   PolyhedralFan(Rays, LS, IncidenceMatrix(Polymake.IncidenceMatrix(Incidence)))
end
function PolyhedralFan(Rays::Union{Oscar.MatElem,AbstractMatrix}, Incidence::Array{Bool,2})
   PolyhedralFan(Rays,IncidenceMatrix(Polymake.IncidenceMatrix(Incidence)))
end

"""
   dim(PF::PolyhedralFan)

Returns the dimension of a polyhedral fan.
"""
dim(PF::PolyhedralFan) = pm_fan(PF).FAN_DIM

"""
   n_maximal_cones(PF::PolyhedralFan)

Returns the number of maximal cones in a polyhedral fan `PF`.
"""
n_maximal_cones(PF::PolyhedralFan) = PF.pm_fan.N_MAXIMAL_CONES

"""
   ambient_dim(PF::PolyhedralFan)

Returns the ambient dimension of a polyhedral fan.
"""
ambient_dim(PF::PolyhedralFan) = pm_fan(PF).FAN_AMBIENT_DIM

"""
   pm_fan(PF::PolyhedralFan)

Get the underlying polymake BigObject
"""
pm_fan(PF::PolyhedralFan) = PF.pm_fan


"""
   lineality_space(PF::PolyhedralFan)

Returns the lineality_space of a polyhedral fan
"""
lineality_space(PF::PolyhedralFan) = pm_fan(PF).LINEALITY_SPACE

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
Base.eltype(::Type{PolyhedralFanRayIterator}) = Polymake.VectorAllocated{Polymake.Rational}
Base.length(iter::PolyhedralFanRayIterator) = n_rays(iter.fan)

"""
   rays(PF::PolyhedralFan)

Returns the rays of a polyhedral fan
"""
rays(PF::PolyhedralFan) = PolyhedralFanRayIterator(PF)

"""
   n_rays(PF::PolyhedralFan)

Returns the number of rays of a polyhedral fan
"""
n_rays(PF::PolyhedralFan) = pm_fan(PF).N_RAYS

"""
   rays_as_point_matrix(PF::PolyhedralFan)

Returns the rays of a polyhedral fan as rows of a matrix
"""
rays_as_point_matrix(PF::PolyhedralFan) = pm_fan(PF).RAYS


"""
   maximal_cones(PF::PolyhedralFan)

Returns the maximal cones of a polyhedral fan
"""

#TODO: should the documentation mention maximal_cones_as_incidence_matrix?
#      similarly for cone ray iterators and facet iterators?
@doc Markdown.doc"""
   maximal_cones(PF::PolyhedralFan, as = :cones)

Returns the maximal cones of the polyhedral fan `PF` in the format defined by `as`.
The allowed values for `as` are
* `cones`: Returns for each maximal cone its realization as a cone
"""
function maximal_cones(PF::PolyhedralFan; as::Symbol = :cones)
   if as == :cones
      MaximalConeIterator(PF)
   else
     throw(ArgumentError("Unsupported `as` argument :" * string(as)))
  end
end

struct MaximalConeIterator
    PF::PolyhedralFan
end

function Base.iterate(iter::MaximalConeIterator, index = 1)
    n_max_cones = n_maximal_cones(iter.PF)
    if index > n_max_cones
        return nothing
    end
    current_cone = Cone(Polymake.fan.cone(iter.PF.pm_fan, index - 1))
    return (current_cone, index + 1)
end
Base.length(iter::MaximalConeIterator) = n_maximal_cones(iter.PF)


function maximal_cones_as_incidence_matrix(PF::PolyhedralFan)
   IncidenceMatrix(PF.pm_fan.MAXIMAL_CONES)
end


"""
   issmooth(PF::PolyhedralFan)

Determine whether the fan is smooth
"""
issmooth(PF::PolyhedralFan) = pm_fan(PF).SMOOTH_FAN

"""
   isregular(PF::PolyhedralFan)

Determine whether the fan is regular, i.e. the normal fan of a polytope
"""
isregular(PF::PolyhedralFan) = pm_fan(PF).REGULAR

"""
   iscomplete(PF::PolyhedralFan)

Determine whether the fan is complete
"""
iscomplete(PF::PolyhedralFan) = pm_fan(PF).COMPLETE

"""
   normal_fan(P::Polyhedron)

Returns the normal fan of a polyhedron
"""
function normal_fan(P::Polyhedron)
   pmp = pm_polytope(P)
   pmnf = Polymake.fan.normal_fan(pmp)
   return PolyhedralFan(pmnf)
end

"""
   face_fan(P::Polyhedron)

Returns the face fan of a polyhedron
"""
function face_fan(P::Polyhedron)
   pmp = pm_polytope(P)
   pmff = Polymake.fan.face_fan(pmp)
   return PolyhedralFan(pmff)
end
