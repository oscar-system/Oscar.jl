########################
# 1: Special attributes of toric varieties
########################

@doc Markdown.doc"""
    chow_ring(v::AbstractNormalToricVariety)

Return the Chow ring of the simplicial and complete toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> ngens(chow_ring(p2))
3
```
"""
@attr MPolyQuo function chow_ring(v::AbstractNormalToricVariety)
    if !issimplicial(v) || !iscomplete(v)
        throw(ArgumentError("The Chow ring is (currently) only supported for simplicial and complete toric varieties."))
    end
    return cohomology_ring(v)
end
export chow_ring
