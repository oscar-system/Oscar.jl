###############################################################################
###############################################################################
### Constructors
###############################################################################
###############################################################################

@doc Markdown.doc"""
    PolyhedralFan{T}(Rays, Cones) where T<:scalar_types

# Arguments
- `R::Matrix`: Rays generating the cones of the fan; encoded row-wise as representative vectors.
- `Cones::IncidenceMatrix`: An incidence matrix; there is a 1 at position (i,j) if cone i has ray j as extremal ray, and 0 otherwise.

A polyhedral fan formed from rays and cones made of these rays. The cones are
given as an IncidenceMatrix, where the columns represent the rays and the rows
represent the cones.

# Examples
To obtain the upper half-space of the plane:
```jldoctest
julia> R = [1 0; 1 1; 0 1; -1 0; 0 -1];

julia> IM=IncidenceMatrix([[1,2],[2,3],[3,4],[4,5],[1,5]]);

julia> PF=PolyhedralFan(R,IM)
A polyhedral fan in ambient dimension 2 with fmpq type coefficients
```
"""
function PolyhedralFan{T}(Rays::Union{SubObjectIterator{<:RayVector}, Oscar.MatElem,AbstractMatrix}, Incidence::IncidenceMatrix) where T<:scalar_types
   PolyhedralFan{T}(Polymake.fan.PolyhedralFan{scalar_type_to_polymake[T]}(
      INPUT_RAYS = Rays,
      INPUT_CONES = Incidence,
   ))
end
function PolyhedralFan{T}(Rays::Union{SubObjectIterator{<:RayVector}, Oscar.MatElem,AbstractMatrix}, LS::Union{SubObjectIterator{<:RayVector}, Oscar.MatElem,AbstractMatrix}, Incidence::IncidenceMatrix) where T<:scalar_types
   PolyhedralFan{T}(Polymake.fan.PolyhedralFan{scalar_type_to_polymake[T]}(
      INPUT_RAYS = Rays,
      INPUT_LINEALITY = LS,
      INPUT_CONES = Incidence,
   ))
end

"""
    pm_object(PF::PolyhedralFan)

Get the underlying polymake object, which can be used via Polymake.jl.
"""
pm_object(PF::PolyhedralFan) = PF.pm_fan

PolyhedralFan(itr::AbstractVector{Cone{T}}) where T<:scalar_types = PolyhedralFan{T}(Polymake.fan.check_fan_objects(pm_object.(itr)...))

#Same construction for when the user gives Matrix{Bool} as incidence matrix
function PolyhedralFan{T}(Rays::Union{SubObjectIterator{<:RayVector}, Oscar.MatElem, AbstractMatrix}, LS::Union{Oscar.MatElem, AbstractMatrix}, Incidence::Matrix{Bool}) where T<:scalar_types
   PolyhedralFan(Rays, LS, IncidenceMatrix(Polymake.IncidenceMatrix(Incidence)))
end
function PolyhedralFan{T}(Rays::Union{SubObjectIterator{<:RayVector}, Oscar.MatElem, AbstractMatrix}, Incidence::Matrix{Bool}) where T<:scalar_types
   PolyhedralFan(Rays, IncidenceMatrix(Polymake.IncidenceMatrix(Incidence)))
end


function PolyhedralFan(C::Cone{T}) where T<:scalar_types
    pmfan = Polymake.fan.check_fan_objects(pm_object(C))
    return PolyhedralFan{T}(pmfan)
end

###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, PF::PolyhedralFan{T}) where T<:scalar_types
    print(io, "A polyhedral fan in ambient dimension $(ambient_dim(PF)) with $T type coefficients")
end
