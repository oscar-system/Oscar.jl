################################################
# 1: Special attributes of toric varieties
################################################

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
    if !is_simplicial(v) || !is_complete(v)
        throw(ArgumentError("The Chow ring is (currently) only supported for simplicial and complete toric varieties"))
    end
    R, _ = PolynomialRing(coefficient_ring(v), coordinate_names(v), cached = false)
    linear_relations = ideal_of_linear_relations(R, v)
    stanley_reisner = stanley_reisner_ideal(R, v)
    return quo(R, linear_relations + stanley_reisner)[1]
end
export chow_ring


@doc Markdown.doc"""
    gens_of_rational_equivalence_classes(v::AbstractNormalToricVariety)

Returns a list of generators of the Chow ring of a complete,
simplicial toric variety.

Recall that the cones of a complete, simplicial toric variety
can be seen as generators of the Chow ring (lemma 12.5.1 in
[CLS11](@cite)). This function first maps each cone to an
element of the Chow ring and then removes elements by taking
rational equivalence into account.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> length(gens_of_rational_equivalence_classes(p2))
6
```
"""
@attr Vector{MPolyQuoElem{fmpq_mpoly}} function gens_of_rational_equivalence_classes(v::AbstractNormalToricVariety)
    g = gens(chow_ring(v))
    r_list = [rays(c) for c in cones(v)]
    return [prod([g[findfirst(x->x==r, rays(v))] for r in rs]) for rs in r_list]
end
export gens_of_rational_equivalence_classes


@doc Markdown.doc"""
    map_gens_of_chow_ring_to_cox_ring(v::AbstractNormalToricVariety)

Returns a dictionary which maps the generators of the chow
ring to monomials in the Cox ring. This dictionary involves
a choice, i.e. is not unique.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> length(map_gens_of_chow_ring_to_cox_ring(p2))
4
```
"""
@attr Dict{fmpq_mpoly, MPolyElem_dec{fmpq, fmpq_mpoly}} function map_gens_of_chow_ring_to_cox_ring(v::AbstractNormalToricVariety)
    g = gens(chow_ring(v))
    g2 = gens(cox_ring(v))
    r_list = [rays(c) for c in cones(v)]
    mapping = Dict{fmpq_mpoly, MPolyElem_dec{fmpq, fmpq_mpoly}}()
    for rs in r_list
      p1 = prod([g[findfirst(x->x==r, rays(v))] for r in rs]).f
      p2 = prod([g2[findfirst(x->x==r, rays(v))] for r in rs])
      coeff = [c for c in AbstractAlgebra.coefficients(p1)][1]
      if coeff != 1
        p1 = 1//coeff * p1
        p2 = 1//coeff * p2
      end
      mapping[p1] = p2
    end
    return mapping
end
export map_gens_of_chow_ring_to_cox_ring
