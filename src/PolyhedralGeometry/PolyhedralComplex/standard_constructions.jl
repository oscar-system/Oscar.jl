# Workaround for turning a PolyhedralFan of polymake into a proper PolyhedralComplex
function polyhedral_complex_workaround(pm::Polymake.BigObject)
    pc = pm
    typename = Polymake.type_name(pm)
    if typename[1:13] == "PolyhedralFan"
        pc = Polymake.fan.PolyhedralComplex(pm)
    end
    typename = Polymake.type_name(pc)
    if typename[1:17] != "PolyhedralComplex"
        error("Input object is not of type PolyhedralFan or PolyhedralComplex")
    end
    fv = Polymake.to_one_based_indexing(pc.FAR_VERTICES)
    mc = pc.MAXIMAL_POLYTOPES
    feasibles = [Polymake.to_zero_based_indexing(Polymake.row(mc, i)) for i in 1:Polymake.nrows(mc) if Polymake.incl(Polymake.row(mc, i), fv)>0]
    return Polymake.fan.PolyhedralComplex(POINTS=pc.VERTICES, INPUT_LINEALITY=pc.LINEALITY_SPACE, INPUT_POLYTOPES=feasibles)
end

@doc Markdown.doc"""
    common_refinement(PC::PolyhedralComplex,PC::PolyhedralComplex)

Return the common refinement of two polyhedral complexes. 

# Example 
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
function common_refinement(PC1::PolyhedralComplex,PC2::PolyhedralComplex) 
    pm_PC1 = pm_object(PC1)
    pm_PC2 = pm_object(PC2)
    result = Polymake.fan.PolyhedralComplex(Polymake.fan.common_refinement(pm_PC1,pm_PC2))
    result = polyhedral_complex_workaround(result)
    return PolyhedralComplex(result)
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
function k_skeleton(PC::PolyhedralComplex,k::Int)
    pm_PC = pm_object(PC)
    ksk = Polymake.fan.PolyhedralComplex(Polymake.fan.k_skeleton(pm_PC,k+1))
    result = polyhedral_complex_workaround(ksk)
    return PolyhedralComplex(result)
end


