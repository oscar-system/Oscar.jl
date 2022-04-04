############################
# Attributes
############################

@doc Markdown.doc"""
    toric_variety(tvs::ToricVanishingSet)

Return the toric variety of the vanishing set `tvs`.
"""
function toric_variety(tvs::ToricVanishingSet)
    return tvs.toric_variety::AbstractNormalToricVarietyss
end
export toric_variety


@doc Markdown.doc"""
    polyhedra(tvs::ToricVanishingSet)

Return the vector of the polyhedra whose complement defines the vanishing set `tvs`.
"""
function polyhedra(tvs::ToricVanishingSet)
    return tvs.ps::Vector{Polyhedron{fmpq}}
end
export polyhedra


@doc Markdown.doc"""
    cohomology_index(tvs::ToricVanishingSet)

Return the cohomology index of the toric vanishing set `tvs`.
"""
function cohomology_index(tvs::ToricVanishingSet)
    return tvs.i::Int
end
export cohomology_index
