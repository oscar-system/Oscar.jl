###############################################################################
###############################################################################
### Definition and constructors
###############################################################################
###############################################################################

struct PolyhedralFan{T} <: PolyhedralObject{T}
  pm_fan::Polymake.BigObject
  parent_field::Field

  PolyhedralFan{T}(pm::Polymake.BigObject, f::Field) where {T<:scalar_types} = new{T}(pm, f)
  PolyhedralFan{QQFieldElem}(pm::Polymake.BigObject) = new{QQFieldElem}(pm, QQ)
end

const _FanLikeType = Union{NormalToricVarietyType,PolyhedralFan}
const _FanLikeTypeQQ = Union{NormalToricVarietyType,PolyhedralFan{QQFieldElem}}

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

julia> IM = incidence_matrix([[1,2],[2,3],[3,4],[4,5],[1,5]]);

julia> PF=polyhedral_fan(IM, R)
Polyhedral fan in ambient dimension 2
```

Polyhedral fan with lineality space:
```jldoctest
julia> R = [1 0 0; 0 0 1];

julia> L = [0 1 0];

julia> IM = incidence_matrix([[1],[2]]);

julia> PF=polyhedral_fan(IM, R, L)
Polyhedral fan in ambient dimension 3

julia> lineality_dim(PF)
1
```
"""
function polyhedral_fan(
  f::scalar_type_or_field,
  Incidence::IncidenceMatrix,
  Rays::AbstractCollection[RayVector],
  LS::Union{AbstractCollection[RayVector],Nothing};
  non_redundant::Bool=false,
)
  parent_field, scalar_type = _determine_parent_and_scalar(f, Rays, LS)
  RM = unhomogenized_matrix(Rays)
  if isnothing(LS)
    LM = Polymake.Matrix{_scalar_type_to_polymake(scalar_type)}(undef, 0, size(RM, 2))
  else
    LM = unhomogenized_matrix(LS)
  end
  if non_redundant
    return PolyhedralFan{scalar_type}(
      Polymake.fan.PolyhedralFan{_scalar_type_to_polymake(scalar_type)}(;
        RAYS=RM, LINEALITY_SPACE=LM, MAXIMAL_CONES=Incidence
      ),
      parent_field,
    )
  else
    return PolyhedralFan{scalar_type}(
      Polymake.fan.PolyhedralFan{_scalar_type_to_polymake(scalar_type)}(;
        INPUT_RAYS=RM, INPUT_LINEALITY=LM, INPUT_CONES=Incidence
      ),
      parent_field,
    )
  end
end
polyhedral_fan(
  f::scalar_type_or_field,
  Incidence::IncidenceMatrix,
  Rays::AbstractCollection[RayVector];
  non_redundant::Bool=false,
) = polyhedral_fan(f, Incidence, Rays, nothing; non_redundant)
polyhedral_fan(
  Incidence::IncidenceMatrix,
  Rays::AbstractCollection[RayVector],
  LS::Union{AbstractCollection[RayVector],Nothing};
  non_redundant::Bool=false,
) = polyhedral_fan(_guess_fieldelem_type(Rays, LS), Incidence, Rays, LS; non_redundant)
polyhedral_fan(
  Incidence::IncidenceMatrix, Rays::AbstractCollection[RayVector]; non_redundant::Bool=false
) = polyhedral_fan(_guess_fieldelem_type(Rays), Incidence, Rays; non_redundant)

"""
    pm_object(PF::PolyhedralFan)

Get the underlying polymake object, which can be used via Polymake.jl.
"""
pm_object(PF::PolyhedralFan) = PF.pm_fan

@doc raw"""
    polyhedral_fan(cones::AbstractVector{Cone{T}}) where T<:scalar_types

Assemble a polyhedral fan from a non-empty list of cones.
"""
function polyhedral_fan(
  cones::AbstractVector{Cone{T}}; non_redundant::Bool=false
) where {T<:scalar_types}
  @req length(cones) > 0 "list of cones must be non-empty"
  if non_redundant
    pmfan = Polymake.fan.fan_from_cones(pm_object.(cones)...)
  else
    pmfan = Polymake.fan.check_fan_objects(pm_object.(cones)...)
  end
  return PolyhedralFan{T}(pmfan, coefficient_field(iterate(cones)[1]))
end

#Same construction for when the user gives Matrix{Bool} as incidence matrix
polyhedral_fan(
  f::scalar_type_or_field,
  Incidence::Matrix{Bool},
  Rays::AbstractCollection[RayVector],
  LS::AbstractCollection[RayVector],
) = polyhedral_fan(f, IncidenceMatrix(Polymake.IncidenceMatrix(Incidence)), Rays, LS)
polyhedral_fan(
  f::scalar_type_or_field, Incidence::Matrix{Bool}, Rays::AbstractCollection[RayVector]
) = polyhedral_fan(f, IncidenceMatrix(Polymake.IncidenceMatrix(Incidence)), Rays)

polyhedral_fan(C::Cone{T}) where {T<:scalar_types} = polyhedral_fan([C])

# shortcut for maximal cones of an existing fan
function polyhedral_fan(iter::SubObjectIterator{Cone{T}}) where {T<:scalar_types}
  if iter.Acc == _maximal_cone && iter.Obj isa PolyhedralFan
    return deepcopy(iter.Obj)
  end
  return polyhedral_fan(collect(iter))
end

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
function polyhedral_fan_from_rays_action(
  f::scalar_type_or_field,
  Rays::AbstractCollection[RayVector],
  MC_reps::IncidenceMatrix,
  perms::AbstractVector{PermGroupElem},
)
  parent_field, scalar_type = _determine_parent_and_scalar(f, Rays)
  pf = Polymake.fan.PolyhedralFan{_scalar_type_to_polymake(scalar_type)}()
  Polymake.take(
    pf,
    "RAYS",
    Polymake.Matrix{_scalar_type_to_polymake(scalar_type)}(unhomogenized_matrix(Rays)),
  )
  d = length(Rays)
  gp = _group_generators_to_pm_arr_arr(perms, d)
  Polymake.take(pf, "GROUP.RAYS_ACTION.GENERATORS", gp)
  Polymake.take(pf, "GROUP.MAXIMAL_CONES_ACTION.MAXIMAL_CONES_GENERATORS", MC_reps)
  return PolyhedralFan{scalar_type}(pf, parent_field)
end
polyhedral_fan_from_rays_action(
  Rays::AbstractCollection[RayVector],
  MC_reps::IncidenceMatrix,
  perms::AbstractVector{PermGroupElem},
) = polyhedral_fan_from_rays_action(QQFieldElem, Rays, MC_reps, perms)

###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, PF::PolyhedralFan{T}) where {T<:scalar_types}
  print(io, "Polyhedral fan in ambient dimension $(ambient_dim(PF))")
  T != QQFieldElem && print(io, " with $T type coefficients")
end
