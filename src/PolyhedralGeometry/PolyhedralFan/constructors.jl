###############################################################################
###############################################################################
### Definition and constructors
###############################################################################
###############################################################################


# We introduce this abstract (hidden) type to allow for other objects to be
# used like polyhedral fans without duplicating too much code, concretely we
# want to be able to directly access rays, maximal_cones, etc for
# NormalToricVariety's.
abstract type _FanLikeType{T} end

struct PolyhedralFan{T} <:_FanLikeType{T}
   pm_fan::Polymake.BigObject
   
   PolyhedralFan{T}(pm::Polymake.BigObject) where T<:scalar_types = new{T}(pm)
end


# Automatic detection of corresponding OSCAR scalar type;
# Avoid, if possible, to increase type stability
PolyhedralFan(p::Polymake.BigObject) = PolyhedralFan{detect_scalar_type(PolyhedralFan, p)}(p)

@doc raw"""
    PolyhedralFan{T}(Rays::AbstractCollection[RayVector], LS::Union{AbstractCollection[RayVector], Nothing}, Incidence::IncidenceMatrix) where T<:scalar_types

Assemble a polyhedral fan from ray generators, lineality generators, and an
`IncidenceMatrix` indicating which rays form a cone.

# Arguments
- `Rays::AbstractCollection[RayVector]`: Rays generating the cones of the fan;
  encoded row-wise as representative vectors.
- `LS::AbstractCollection[RayVector]`: Contains row-wise generators of the
  lineality space of the fan. (optional argument)
- `Cones::IncidenceMatrix`: An incidence matrix; there is a 1 at position (i,j)
  if cone i has ray j as extremal ray, and 0 otherwise.


# Examples
To obtain the upper half-space of the plane:
```jldoctest
julia> R = [1 0; 1 1; 0 1; -1 0; 0 -1];

julia> IM=IncidenceMatrix([[1,2],[2,3],[3,4],[4,5],[1,5]]);

julia> PF=PolyhedralFan(R,IM)
Polyhedral fan in ambient dimension 2
```

Polyhedral fan with lineality space:
```jldoctest
julia> R = [1 0 0; 0 0 1];

julia> L = [0 1 0];

julia> IM = IncidenceMatrix([[1],[2]]);

julia> PF=PolyhedralFan(R, L, IM)
Polyhedral fan in ambient dimension 3

julia> lineality_dim(PF)
1
```
"""
function PolyhedralFan{T}(Rays::AbstractCollection[RayVector], LS::Union{AbstractCollection[RayVector], Nothing}, Incidence::IncidenceMatrix; non_redundant::Bool = false) where T<:scalar_types
  RM = unhomogenized_matrix(Rays)
  if isnothing(LS)
    LM = Polymake.Matrix{scalar_type_to_polymake[T]}(undef, 0, size(RM, 2))
  else
    LM = unhomogenized_matrix(LS)
  end
  if non_redundant
    return PolyhedralFan{T}(Polymake.fan.PolyhedralFan{scalar_type_to_polymake[T]}(
      RAYS = RM,
      LINEALITY_SPACE = LM,
      MAXIMAL_CONES = Incidence,
    ))
  else
    return PolyhedralFan{T}(Polymake.fan.PolyhedralFan{scalar_type_to_polymake[T]}(
      INPUT_RAYS = RM,
      INPUT_LINEALITY = LM,
      INPUT_CONES = Incidence,
    ))
  end
end
PolyhedralFan{T}(Rays::AbstractCollection[RayVector], Incidence::IncidenceMatrix; non_redundant::Bool = false) where T<:scalar_types = PolyhedralFan{T}(Rays, nothing, Incidence; non_redundant=non_redundant)
PolyhedralFan(Rays::AbstractCollection[RayVector], LS::Union{AbstractCollection[RayVector], Nothing}, Incidence::IncidenceMatrix; non_redundant::Bool = false) = PolyhedralFan{QQFieldElem}(Rays, LS, Incidence; non_redundant = non_redundant)
PolyhedralFan(Rays::AbstractCollection[RayVector], Incidence::IncidenceMatrix; non_redundant::Bool = false) = PolyhedralFan{QQFieldElem}(Rays, Incidence; non_redundant = non_redundant)

"""
    pm_object(PF::PolyhedralFan)

Get the underlying polymake object, which can be used via Polymake.jl.
"""
pm_object(PF::PolyhedralFan) = PF.pm_fan

PolyhedralFan(itr::AbstractVector{Cone{T}}) where T<:scalar_types = PolyhedralFan{T}(Polymake.fan.check_fan_objects(pm_object.(itr)...))

#Same construction for when the user gives Matrix{Bool} as incidence matrix
function PolyhedralFan{T}(Rays::AbstractCollection[RayVector], LS::AbstractCollection[RayVector], Incidence::Matrix{Bool}) where T<:scalar_types
   PolyhedralFan(Rays, LS, IncidenceMatrix(Polymake.IncidenceMatrix(Incidence)))
end
function PolyhedralFan{T}(Rays::AbstractCollection[RayVector], Incidence::Matrix{Bool}) where T<:scalar_types
   PolyhedralFan(Rays, IncidenceMatrix(Polymake.IncidenceMatrix(Incidence)))
end


function PolyhedralFan(C::Cone{T}) where T<:scalar_types
    pmfan = Polymake.fan.check_fan_objects(pm_object(C))
    return PolyhedralFan{T}(pmfan)
end

###############################################################################
###############################################################################
### From group action on maximal cones
###############################################################################
###############################################################################
@doc raw"""
    polyhedral_fan_from_rays_action(::Type{T}, Rays::AbstractCollection[RayVector], MC_reps::IncidenceMatrix, perms::AbstractVector{PermGroupElem}) where T<:scalar_types

Construct a polyhedral fan with a group action.

# Arguments
- `Rays`: The rays of the fan
- `MC_reps`: `IncidenceMatrix` whose rows give the indices of the rays forming
  representatives of the maximal cones under the group action.
- `perms`: A vector of permutations `PermGroupElem` that form generators of the
  group acting on the rays of the fan.

"""
function polyhedral_fan_from_rays_action(::Type{T}, Rays::AbstractCollection[RayVector], MC_reps::IncidenceMatrix, perms::AbstractVector{PermGroupElem}) where T<:scalar_types
    pf = Polymake.fan.PolyhedralFan()
    Polymake.take(pf, "RAYS", Polymake.Matrix(unhomogenized_matrix(Rays)))
    d = length(Rays)
    gp = _group_generators_to_pm_arr_arr(perms, d)
    Polymake.take(pf, "GROUP.REPRESENTATIVE_MAXIMAL_CONES", MC_reps)
    Polymake.take(pf, "GROUP.RAYS_ACTION.GENERATORS", gp)
    return PolyhedralFan{T}(pf)
end
polyhedral_fan_from_rays_action(Rays::AbstractCollection[RayVector], MC_reps::IncidenceMatrix, perms::AbstractVector{PermGroupElem}) = polyhedral_fan_from_rays_action(QQFieldElem, Rays, MC_reps, perms)


###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, PF::PolyhedralFan{T}) where T<:scalar_types
    print(io, "Polyhedral fan in ambient dimension $(ambient_dim(PF))")
    T != QQFieldElem && print(io, " with $T type coefficients")
end
