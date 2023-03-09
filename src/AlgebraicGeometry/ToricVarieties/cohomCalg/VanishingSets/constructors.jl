#########################
# 1: The Julia type for toric vanishing sets
#########################

@attributes mutable struct ToricVanishingSet
    toric_variety::AbstractNormalToricVariety
    ps::Vector{Polyhedron{QQFieldElem}}
    i::Int
    function ToricVanishingSet(toric_variety::AbstractNormalToricVariety, ps::Vector{Polyhedron{QQFieldElem}}, i::Int)
        if !all(p -> ambient_dim(p) == rank(picard_group(toric_variety)), ps)
            throw(ArgumentError("The ambient dimensions of the polyhedra must match the rank as the picard group of the toric variety"))
        end
        if !(0 <= i <= dim(toric_variety))
            throw(ArgumentError("The cohomology index must not be negative and not larger than $(dim(toric_variety))"))
        end
        return new(toric_variety, ps, i)
    end
end

toric_vanishing_set(v::AbstractNormalToricVariety, ps::Vector{Polyhedron{QQFieldElem}}, i::Int) = ToricVanishingSet(v, ps)


######################
# 2: Display
######################s

function Base.show(io::IO, tvs::ToricVanishingSet)
    join(io, "Toric vanishing set for cohomology index $(cohomology_index(tvs))")
end
