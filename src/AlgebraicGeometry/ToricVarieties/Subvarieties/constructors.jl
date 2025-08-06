##############################################################
# 1: The Julia type for subvarieties of normal toric varieties
##############################################################

@attributes mutable struct ClosedSubvarietyOfToricVariety
    toric_variety::NormalToricVarietyType
    defining_ideal::MPolyIdeal
    function ClosedSubvarietyOfToricVariety(toric_variety::NormalToricVarietyType, defining_ideal::MPolyIdeal)
      @req is_simplicial(toric_variety) "Currently, closed subvarieties are only supported for simplicial toric varieties"
      @req base_ring(defining_ideal) == cox_ring(toric_variety) "The defining ideal must be contained in the Cox ring of the toric supervariety"
      return new(toric_variety, defining_ideal)
    end
end


##############################################################
# 2: Generic constructor
##############################################################

@doc raw"""
    closed_subvariety_of_toric_variety(toric_variety::NormalToricVarietyType, defining_polynomials::Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}})

Construct the closed subvariety of a simplicial normal toric variety.
The defining data for the closed subvariety is a list of homogeneous
polynomials, all of which must be elements of the Cox ring of the
toric variety in question. The common vanishing locus of these polynomials
defines the closed subvariety in question. By proposition 5.2.4 in
[CLS11](@cite) every closed subvariety of a simplicial toric variety
arises in this way.

# Examples
```jldoctest
julia> f2 = hirzebruch_surface(NormalToricVariety, 2);

julia> (t1, x1, t2, x2) = gens(cox_ring(f2));

julia> closed_subvariety_of_toric_variety(f2, [t1])
Closed subvariety of a normal toric variety
```
"""
closed_subvariety_of_toric_variety(toric_variety::NormalToricVarietyType, defining_polynomials::Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}) = ClosedSubvarietyOfToricVariety(toric_variety, ideal(defining_polynomials))


@doc raw"""
    closed_subvariety_of_toric_variety(toric_variety::NormalToricVarietyType, defining_ideal::MPolyIdeal)

Construct the closed subvariety of a simplicial normal toric variety.
The defining data for the closed subvariety is an ideal of the Cox ring of the
toric variety in question. By proposition 5.2.4 in
[CLS11](@cite) every closed subvariety of a simplicial toric variety
arises in this way.

# Examples
```jldoctest
julia> f2 = hirzebruch_surface(NormalToricVariety, 2);

julia> (t1, x1, t2, x2) = gens(cox_ring(f2));

julia> closed_subvariety_of_toric_variety(f2, ideal([t1]))
Closed subvariety of a normal toric variety
```
"""
closed_subvariety_of_toric_variety(toric_variety::NormalToricVarietyType, defining_ideal::MPolyIdeal) = ClosedSubvarietyOfToricVariety(toric_variety, defining_ideal)


######################
# 3: Display
######################s

function Base.show(io::IO, c::ClosedSubvarietyOfToricVariety)
    properties_string = ["Closed"]
    push_attribute_if_exists!(properties_string, c, :is_empty, "empty")
    push!(properties_string, "subvariety of a normal toric variety")
    join(io, properties_string, ", ", " ")
end
