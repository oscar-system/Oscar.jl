@doc Markdown.doc"""
    common_refinement(PC1::PolyhedralComplex{T},PC2::PolyhedralComplex{T}) where T<:scalar_types

Return the common refinement of two polyhedral complexes. 

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2,3]])
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

julia> PC1 = PolyhedralComplex(IM,VR1)
A polyhedral complex in ambient dimension 2

julia> PC2 = PolyhedralComplex(IM,VR2)
A polyhedral complex in ambient dimension 2

julia> common_refinement(PC1,PC2)
A polyhedral complex in ambient dimension 2
```
"""
function common_refinement(PC1::PolyhedralComplex{T},PC2::PolyhedralComplex{T}) where T<:scalar_types
    pm_PC1 = pm_object(PC1)
    pm_PC2 = pm_object(PC2)
    result = Polymake.fan.PolyhedralComplex{scalar_type_to_polymake[T]}(Polymake.fan.common_refinement(pm_PC1,pm_PC2))
    return PolyhedralComplex{T}(result)
end


@doc Markdown.doc"""
     k_skeleton(PC::PolyhedralComplex,k::Int)

Return the k-skeleton of a polyhedral complex.
#Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2,3]])
1×3 IncidenceMatrix
[1, 2, 3]

julia> VR = [0 0; 1 0; 1 1]
3×2 Matrix{Int64}:
 0  0
 1  0
 1  1

julia> PC1 = PolyhedralComplex(IM,VR)
A polyhedral complex in ambient dimension 2

julia> k_skeleton(PC1,1)
A polyhedral complex in ambient dimension 2
```
"""
function k_skeleton(PC::PolyhedralComplex{T},k::Int) where T<:scalar_types
    pm_PC = pm_object(PC)
    ksk = Polymake.fan.PolyhedralComplex{scalar_type_to_polymake[T]}(Polymake.fan.k_skeleton(pm_PC,k))
    return PolyhedralComplex{T}(ksk)
end


