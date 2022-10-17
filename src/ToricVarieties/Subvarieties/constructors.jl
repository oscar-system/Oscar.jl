##############################################################
# 1: The Julia type for subvarieties of normal toric varieties
##############################################################

@attributes mutable struct ClosedSubvarietyOfToricVariety
    toric_variety::AbstractNormalToricVariety
    defining_ideal::MPolyIdeal
    function ClosedSubvarietyOfToricVariety(toric_variety::AbstractNormalToricVariety, defining_ideal::MPolyIdeal)
      if !is_simplicial(toric_variety)
        throw(ArgumentError("Currently, closed subvarieties are only supported for simplicial toric varieties"))
      end
      if parent(gens(defining_ideal)[1]) != cox_ring(toric_variety)
        throw(ArgumentError("The defining ideal must be contained in the Cox ring of the toric supervariety"))
      end
      return new(toric_variety, defining_ideal)
    end
end
export ClosedSubvarietyOfToricVariety


##############################################################
# 2: Generic constructor
##############################################################

@doc Markdown.doc"""
    ClosedSubvarietyOfToricVariety(toric_variety::AbstractNormalToricVariety, defining_polynomials::Vector{MPolyElem_dec{fmpq, fmpq_mpoly}})

Construct the closed subvariety of a simplicial normal toric variety.
The defining data for the closed subvariety is a list of homogeneous
polynomials, all of which must be elements of the Cox ring of the
toric variety in question. The common vanishing locus of these polynomials
defines the closed subvariety in question. By proposition 5.2.4 in
[CLS11](@cite) every closed subvariety of a simplicial toric variety
arises in this way.

# Examples
```jldoctest
julia> f2 = hirzebruch_surface(2);

julia> (t1, x1, t2, x2) = gens(cox_ring(f2));

julia> ClosedSubvarietyOfToricVariety(f2, [t1])
A closed subvariety of a normal toric variety
```
"""
ClosedSubvarietyOfToricVariety(toric_variety::AbstractNormalToricVariety, defining_polynomials::Vector{MPolyElem_dec{fmpq, fmpq_mpoly}}) = ClosedSubvarietyOfToricVariety(toric_variety, ideal(defining_polynomials))
export ClosedSubvarietyOfToricVariety


######################
# 3: Display
######################s

function Base.show(io::IO, c::ClosedSubvarietyOfToricVariety)
    properties_string = ["A closed"]
    push_attribute_if_exists!(properties_string, c, :is_empty, "empty")
    push!(properties_string, "subvariety of a normal toric variety")
    join(io, properties_string, ", ", " ")
end
