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
## Scalar product
###############################################################################

function *(c::QQFieldElem, Sigma::PolyhedralComplex)
    # if scalar is zero, return polyhedral complex consisting only of the origin
    if iszero(c)
        return polyhedral_complex(convex_hull(zero_matrix(QQ,1,ambient_dim(Sigma))))
    end

    # if scalar is non-zero, multiple all vertices and rays by said scalar
    SigmaVertsAndRays = vertices_and_rays(Sigma)
    SigmaRayIndices = findall(vr -> vr isa RayVector, SigmaVertsAndRays)
    SigmaLineality = lineality_space(Sigma)
    SigmaIncidence = maximal_polyhedra(IncidenceMatrix,Sigma)
    return polyhedral_complex(SigmaIncidence, multiply_by_nonzero_scalar.(SigmaVertsAndRays,c), SigmaRayIndices, SigmaLineality)
end
*(c::RationalUnion, Sigma::PolyhedralComplex) = QQ(c)*Sigma
*(Sigma::PolyhedralComplex, c::RationalUnion) = c*Sigma


###############################################################################
## Negation
###############################################################################

function -(Sigma::PolyhedralComplex)
    SigmaVertsAndRays = vertices_and_rays(Sigma)
    SigmaRayIndices = findall(vr -> vr isa RayVector, SigmaVertsAndRays)
    SigmaLineality = lineality_space(Sigma)
    SigmaIncidence = maximal_polyhedra(IncidenceMatrix,Sigma)
    return polyhedral_complex(SigmaIncidence, -SigmaVertsAndRays, SigmaRayIndices, SigmaLineality)
end

###############################################################################
## Translation
###############################################################################

function translate_by_vector(u::PointVector{QQFieldElem}, v::Vector{QQFieldElem})
    return u .+ v
end
function translate_by_vector(u::RayVector{QQFieldElem}, ::Vector{QQFieldElem})
    return u
end
function +(v::Vector{QQFieldElem}, Sigma::PolyhedralComplex)
    @req length(v)==ambient_dim(Sigma) "ambient dimension mismatch"
    SigmaVertsAndRays = vertices_and_rays(Sigma)
    SigmaRayIndices = findall(vr -> vr isa RayVector, SigmaVertsAndRays)
    SigmaLineality = lineality_space(Sigma)
    SigmaIncidence = maximal_polyhedra(IncidenceMatrix,Sigma)
    return polyhedral_complex(SigmaIncidence, translate_by_vector.(SigmaVertsAndRays,Ref(v)), SigmaRayIndices, SigmaLineality)
end
+(v::Vector{ZZRingElem}, Sigma::PolyhedralComplex) = QQ.(v)+Sigma
+(v::Vector{Rational}, Sigma::PolyhedralComplex) = QQ.(v)+Sigma
+(v::Vector{Int}, Sigma::PolyhedralComplex) = QQ.(v)+Sigma

+(Sigma::PolyhedralComplex, v::Vector{QQFieldElem}) = v+Sigma
+(Sigma::PolyhedralComplex, v::Vector{<:RationalUnion}) = QQ.(v)+Sigma

# Vector addition for polyhedral fans
+(Sigma::PolyhedralFan, v::Vector) = polyhedral_complex(Sigma)+v
+(v::Vector, Sigma::PolyhedralFan) = v+polyhedral_complex(Sigma)
