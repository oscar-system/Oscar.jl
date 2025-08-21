#########################
# 1: The Julia type for toric vanishing sets
#########################

@attributes mutable struct ToricVanishingSet
    toric_variety::NormalToricVarietyType
    ps::Vector{Polyhedron{QQFieldElem}}
    cis::Vector{Int}
    function ToricVanishingSet(toric_variety::NormalToricVarietyType, ps::Vector{Polyhedron{QQFieldElem}}, cis::Vector{Int})
        if !all(p -> ambient_dim(p) == torsion_free_rank(picard_group_with_map(toric_variety)[1]), ps)
            throw(ArgumentError("The ambient dimensions of the polyhedra must match the rank as the picard group of the toric variety"))
        end
        if !all(i -> 0 <= i <= dim(toric_variety), cis)
            throw(ArgumentError("The cohomology indices must not be negative and not larger than $(dim(toric_variety))"))
        end
        return new(toric_variety, ps, cis)
    end
end

toric_vanishing_set(v::NormalToricVarietyType, ps::Vector{Polyhedron{QQFieldElem}}, cis::Vector{Int}) = ToricVanishingSet(v, ps, cis)


######################
# 2: Display
######################s

function Base.show(io::IO, tvs::ToricVanishingSet)
    join(io, "Toric vanishing set for cohomology indices $(cohomology_indices(tvs))")
end
