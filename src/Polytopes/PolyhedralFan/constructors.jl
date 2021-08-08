###############################################################################
###############################################################################
### Definition and constructors
###############################################################################
###############################################################################

@doc Markdown.doc"""
    PolyhedralFan(Rays, Cones)

A polyhedral fan formed from rays and cones made of these rays. The cones are
given as an IncidenceMatrix, where the columns represent the rays and the rows
represent the cones. There is a 1 at position (i,j) if cone i has ray j as
extremal ray, otherwise there is a 0.
"""
struct PolyhedralFan
   pm_fan::Polymake.BigObject
   function PolyhedralFan(pm::Polymake.BigObject)
      return new(pm)
   end
end
function PolyhedralFan(Rays::Union{Oscar.MatElem,AbstractMatrix}, Incidence::IncidenceMatrix)
   arr = @Polymake.convert_to Array{Set{Int}} Polymake.common.rows(Incidence.pm_incidencematrix)
   PolyhedralFan(Polymake.fan.PolyhedralFan{Polymake.Rational}(
      INPUT_RAYS = matrix_for_polymake(Rays),
      INPUT_CONES = arr,
   ))
end
function PolyhedralFan(Rays::Union{Oscar.MatElem,AbstractMatrix}, LS::Union{Oscar.MatElem,AbstractMatrix}, Incidence::IncidenceMatrix)
   arr = @Polymake.convert_to Array{Set{Int}} Polymake.common.rows(Incidence.pm_incidencematrix)
   PolyhedralFan(Polymake.fan.PolyhedralFan{Polymake.Rational}(
      INPUT_RAYS = matrix_for_polymake(Rays),
      INPUT_LINEALITY = matrix_for_polymake(LS),
      INPUT_CONES = arr,
   ))
end

"""
    pm_fan(PF::PolyhedralFan)

Get the underlying polymake object, which can be used via Polymake.jl.
"""
pm_fan(PF::PolyhedralFan) = PF.pm_fan


function PolyhedralFan(itr)
   cones = collect(Cone, itr)
   BigObjectArray = Polymake.Array{Polymake.BigObject}(length(cones))
   for i in 1:length(cones)
      BigObjectArray[i] = pm_cone(cones[i])
   end
   PolyhedralFan(Polymake.fan.check_fan_objects(BigObjectArray))
end



#Same construction for when the user gives Matrix{Bool} as incidence matrix
function PolyhedralFan(Rays::Union{Oscar.MatElem,AbstractMatrix}, LS::Union{Oscar.MatElem,AbstractMatrix}, Incidence::Matrix{Bool})
   PolyhedralFan(Rays, LS, IncidenceMatrix(Polymake.IncidenceMatrix(Incidence)))
end
function PolyhedralFan(Rays::Union{Oscar.MatElem,AbstractMatrix}, Incidence::Matrix{Bool})
   PolyhedralFan(Rays,IncidenceMatrix(Polymake.IncidenceMatrix(Incidence)))
end

###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, PF::PolyhedralFan)
    print(io, "A polyhedral fan in ambient dimension $(ambient_dim(PF))")
end

