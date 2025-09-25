@doc raw"""
    common_refinement(PC1::PolyhedralComplex{T},PC2::PolyhedralComplex{T}) where T<:scalar_types

Return the common refinement of two polyhedral complexes.

# Examples
```jldoctest
julia> IM = incidence_matrix([[1,2,3]])
1×3 IncidenceMatrix
 [1, 2, 3]

julia> VR1 = [0 0; 1 0; 1 1]
3×2 Matrix{Int64}:
 0  0
 1  0
 1  1

julia> VR2 = [0 0; 0 1; 1 1]
3×2 Matrix{Int64}:
 0  0
 0  1
 1  1

julia> PC1 = polyhedral_complex(IM,VR1)
Polyhedral complex in ambient dimension 2

julia> PC2 = polyhedral_complex(IM,VR2)
Polyhedral complex in ambient dimension 2

julia> common_refinement(PC1,PC2)
Polyhedral complex in ambient dimension 2
```
"""
function common_refinement(
  PC1::PolyhedralComplex{T}, PC2::PolyhedralComplex{T}
) where {T<:scalar_types}
  U, f = _promote_scalar_field(coefficient_field(PC1), coefficient_field(PC2))
  pm_PC1 = pm_object(PC1)
  pm_PC2 = pm_object(PC2)
  result = Polymake.fan.PolyhedralComplex{_scalar_type_to_polymake(T)}(
    Polymake.fan.common_refinement(pm_PC1, pm_PC2)
  )
  return PolyhedralComplex{T}(result, f)
end

@doc raw"""
    k_skeleton(PC::PolyhedralComplex,k::Int)

Return the k-skeleton of a polyhedral complex.

# Examples
```jldoctest
julia> IM = incidence_matrix([[1,2,3]])
1×3 IncidenceMatrix
 [1, 2, 3]

julia> VR = [0 0; 1 0; 1 1]
3×2 Matrix{Int64}:
 0  0
 1  0
 1  1

julia> PC1 = polyhedral_complex(IM,VR)
Polyhedral complex in ambient dimension 2

julia> k_skeleton(PC1,1)
Polyhedral complex in ambient dimension 2
```
"""
function k_skeleton(PC::PolyhedralComplex{T}, k::Int) where {T<:scalar_types}
  pm_PC = pm_object(PC)
  ksk = Polymake.fan.PolyhedralComplex{_scalar_type_to_polymake(T)}(
    Polymake.fan.k_skeleton(pm_PC, k)
  )
  return PolyhedralComplex{T}(ksk, coefficient_field(PC))
end

###############################################################################
## Scalar multiplication
###############################################################################

# helper functions to ensure that point and ray vectors scale differently:
# - point vectors represent points in space and should be scaled normally
# - ray vectors represent directions and they should only be scaled by 0 or +/-1
function scale(
  u::PointVector{<:scalar_types_extended}, c::scalar_types_extended
)
  return u * c
end
function scale(
  u::RayVector{<:scalar_types_extended}, c::scalar_types_extended
)
  if is_positive(c)
    return u
  end
  if is_negative(c)
    return -u
  end
  return 0*u
end

_empty_complex(cf, dim) = polyhedral_complex(cf, incidence_matrix(0, 0), zero_matrix(cf, 0, dim); non_redundant=true)
_origin_complex(cf, dim) = polyhedral_complex(cf, incidence_matrix(1, 1, [[1]]), zero_matrix(cf, 1, dim); non_redundant=true)

function *(c::scalar_types_extended, Sigma::PolyhedralComplex)
  n_maximal_polyhedra(Sigma) == 0 && return _empty_complex(coefficient_field(Sigma), ambient_dim(Sigma))
  iszero(c) && return _origin_complex(coefficient_field(Sigma), ambient_dim(Sigma))

  # if scalar is non-zero, multiple all vertices and rays by said scalar
  SigmaVertsAndRays = vertices_and_rays(Sigma)
  SigmaRayIndices = findall(vr -> vr isa RayVector, SigmaVertsAndRays)
  SigmaLineality = lineality_space(Sigma)
  SigmaIncidence = maximal_polyhedra(IncidenceMatrix, Sigma)
  return polyhedral_complex(
    coefficient_field(Sigma),
    SigmaIncidence,
    scale.(SigmaVertsAndRays,c),
    SigmaRayIndices,
    SigmaLineality;
    non_redundant=true
  )
end
*(Sigma::PolyhedralComplex, c::scalar_types_extended) = c * Sigma

###############################################################################
## Negation
###############################################################################

function -(Sigma::PolyhedralComplex)
  SigmaVertsAndRays = vertices_and_rays(Sigma)
  SigmaRayIndices = findall(vr -> vr isa RayVector, SigmaVertsAndRays)
  SigmaLineality = lineality_space(Sigma)
  SigmaIncidence = maximal_polyhedra(IncidenceMatrix, Sigma)
  return polyhedral_complex(
    coefficient_field(Sigma), SigmaIncidence, -SigmaVertsAndRays, SigmaRayIndices, SigmaLineality; non_redundant=true
  )
end

###############################################################################
## Translation
###############################################################################

# helper functions to ensure that point and ray vectors are translated differently:
# - point vectors represent points in space and should be translated normally
# - ray vectors represent directions and they should not change at all
function translate(
  u::PointVector{<:scalar_types_extended}, v::AbstractVector{<:scalar_types_extended}
)
  return u + v
end
function translate(
  u::RayVector{<:scalar_types_extended}, ::AbstractVector{<:scalar_types_extended}
)
  return u
end
function +(v::Vector{<:scalar_types_extended}, Sigma::PolyhedralComplex)
  @req length(v) == ambient_dim(Sigma) "ambient dimension mismatch"
  SigmaVertsAndRays = vertices_and_rays(Sigma)
  SigmaRayIndices = findall(vr -> vr isa RayVector, SigmaVertsAndRays)
  SigmaLineality = lineality_space(Sigma)
  SigmaIncidence = maximal_polyhedra(IncidenceMatrix, Sigma)
  return polyhedral_complex(
    coefficient_field(Sigma),
    SigmaIncidence,
    translate.(SigmaVertsAndRays, Ref(v)),
    SigmaRayIndices,
    SigmaLineality,
  )
end
+(Sigma::PolyhedralComplex, v::Vector{<:scalar_types_extended}) = v + Sigma

# Vector addition for polyhedral fans
+(Sigma::PolyhedralFan, v::Vector{<:scalar_types_extended}) = polyhedral_complex(Sigma) + v
+(v::Vector{<:scalar_types_extended}, Sigma::PolyhedralFan) = v + polyhedral_complex(Sigma)
