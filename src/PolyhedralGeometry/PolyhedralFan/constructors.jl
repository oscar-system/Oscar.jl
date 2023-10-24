###############################################################################
###############################################################################
### Definition and constructors
###############################################################################
###############################################################################

struct PolyhedralFan{T} <: PolyhedralObject{T}
    pm_fan::Polymake.BigObject
    parent_field::Field

    PolyhedralFan{T}(pm::Polymake.BigObject, f::Field) where T<:scalar_types = new{T}(pm, f)
    PolyhedralFan{QQFieldElem}(pm::Polymake.BigObject) = new{QQFieldElem}(pm, QQ)
end

const _FanLikeType = Union{NormalToricVarietyType, PolyhedralFan}
const _FanLikeTypeQQ = Union{NormalToricVarietyType, PolyhedralFan{QQFieldElem}}

# Automatic detection of corresponding OSCAR scalar type;
# Avoid, if possible, to increase type stability
function polyhedral_fan(p::Polymake.BigObject)
    T, f = _detect_scalar_and_field(PolyhedralFan, p)
    return PolyhedralFan{T}(p, f)
end

@doc raw"""
    polyhedral_fan(T, Rays::AbstractCollection[RayVector], LS::Union{AbstractCollection[RayVector], Nothing}, Incidence::IncidenceMatrix) where T<:scalar_types

Assemble a polyhedral fan from ray generators, lineality generators, and an
`IncidenceMatrix` indicating which rays form a cone.

# Arguments
- `T`: `Type` or parent `Field` of scalar to use, defaults to `QQFieldElem`.
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

julia> PF=polyhedral_fan(R,IM)
Polyhedral fan in ambient dimension 2
```

Polyhedral fan with lineality space:
```jldoctest
julia> R = [1 0 0; 0 0 1];

julia> L = [0 1 0];

julia> IM = IncidenceMatrix([[1],[2]]);

julia> PF=polyhedral_fan(R, L, IM)
Polyhedral fan in ambient dimension 3

julia> lineality_dim(PF)
1
```
"""
function polyhedral_fan(f::scalar_type_or_field,
                        Rays::AbstractCollection[RayVector], 
                        LS::Union{AbstractCollection[RayVector], Nothing},
                        Incidence::IncidenceMatrix; 
                        non_redundant::Bool = false)
  parent_field, scalar_type = _determine_parent_and_scalar(f, Rays, LS)
  RM = unhomogenized_matrix(Rays)
  if isnothing(LS)
    LM = Polymake.Matrix{_scalar_type_to_polymake(scalar_type)}(undef, 0, size(RM, 2))
  else
    LM = unhomogenized_matrix(LS)
  end
  if non_redundant
    return PolyhedralFan{scalar_type}(Polymake.fan.PolyhedralFan{_scalar_type_to_polymake(scalar_type)}(
      RAYS = RM,
      LINEALITY_SPACE = LM,
      MAXIMAL_CONES = Incidence,
    ), parent_field)
  else
    return PolyhedralFan{scalar_type}(Polymake.fan.PolyhedralFan{_scalar_type_to_polymake(scalar_type)}(
      INPUT_RAYS = RM,
      INPUT_LINEALITY = LM,
      INPUT_CONES = Incidence,
    ), parent_field)
  end
end
polyhedral_fan(f::scalar_type_or_field, Rays::AbstractCollection[RayVector], Incidence::IncidenceMatrix; non_redundant::Bool = false) = polyhedral_fan(f, Rays, nothing, Incidence; non_redundant=non_redundant)
polyhedral_fan(Rays::AbstractCollection[RayVector], LS::Union{AbstractCollection[RayVector], Nothing}, Incidence::IncidenceMatrix; non_redundant::Bool = false) = polyhedral_fan(QQFieldElem, Rays, LS, Incidence; non_redundant = non_redundant)
polyhedral_fan(Rays::AbstractCollection[RayVector], Incidence::IncidenceMatrix; non_redundant::Bool = false) = polyhedral_fan(QQFieldElem, Rays, Incidence; non_redundant = non_redundant)

"""
    pm_object(PF::PolyhedralFan)

Get the underlying polymake object, which can be used via Polymake.jl.
"""
pm_object(PF::PolyhedralFan) = PF.pm_fan

function polyhedral_fan(itr::AbstractVector{Cone{T}}) where T<:scalar_types
  @req length(itr) > 0 "list of cones must be non-empty"
  pmfan = Polymake.fan.check_fan_objects(pm_object.(itr)...)
  if pmfan.N_MAXIMAL_CONES == 0
    # if check fan returns no cones then the rays are empty and we have just one trivial cone (maybe with lineality)
    C = itr[1]
    return polyhedral_fan(coefficient_field(C), rays(C), lineality_space(C), IncidenceMatrix(1,0); non_redundant=true)
  end
  return PolyhedralFan{T}(pmfan, coefficient_field(iterate(itr)[1]))
end

#Same construction for when the user gives Matrix{Bool} as incidence matrix
polyhedral_fan(f::scalar_type_or_field, Rays::AbstractCollection[RayVector], LS::AbstractCollection[RayVector], Incidence::Matrix{Bool}) =
  polyhedral_fan(f, Rays, LS, IncidenceMatrix(Polymake.IncidenceMatrix(Incidence)))
polyhedral_fan(f::scalar_type_or_field, Rays::AbstractCollection[RayVector], Incidence::Matrix{Bool}) =
  polyhedral_fan(f, Rays, IncidenceMatrix(Polymake.IncidenceMatrix(Incidence)))

polyhedral_fan(C::Cone{T}) where T<:scalar_types = polyhedral_fan([C])

###############################################################################
###############################################################################
### From group action on maximal cones
###############################################################################
###############################################################################
@doc raw"""
    polyhedral_fan_from_rays_action([::Union{Type{T}, Field} = QQFieldElem,] Rays::AbstractCollection[RayVector], MC_reps::IncidenceMatrix, perms::AbstractVector{PermGroupElem}) where T<:scalar_types

Construct a polyhedral fan with a group action.

# Arguments
- The first argument either specifies the `Type` of its coefficients or their
parent `Field`.
- `Rays`: The rays of the fan
- `MC_reps`: `IncidenceMatrix` whose rows give the indices of the rays forming
  representatives of the maximal cones under the group action.
- `perms`: A vector of permutations `PermGroupElem` that form generators of the
  group acting on the rays of the fan.

"""
function polyhedral_fan_from_rays_action(f::scalar_type_or_field, Rays::AbstractCollection[RayVector], MC_reps::IncidenceMatrix, perms::AbstractVector{PermGroupElem})
  parent_field, scalar_type = _determine_parent_and_scalar(f, Rays)
  pf = Polymake.fan.PolyhedralFan{_scalar_type_to_polymake(scalar_type)}()
  Polymake.take(pf, "RAYS", Polymake.Matrix{_scalar_type_to_polymake(scalar_type)}(unhomogenized_matrix(Rays)))
  d = length(Rays)
  gp = _group_generators_to_pm_arr_arr(perms, d)
  Polymake.take(pf, "GROUP.RAYS_ACTION.GENERATORS", gp)
  Polymake.take(pf, "GROUP.MAXIMAL_CONES_ACTION.MAXIMAL_CONES_GENERATORS", MC_reps)
  return PolyhedralFan{scalar_type}(pf, parent_field)
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
