################################################
# 1: Special attributes of toric varieties
################################################

@doc raw"""
    chow_ring(v::NormalToricVarietyType)

Return the Chow ring of the simplicial toric variety `v`.

While [CLS11](@cite) focus on simplicial and complete varieties to
define the Chow ring, it was described in [Peg14](@cite) that this
notion can also be extended to non-complete varieties. We explicitly
support the Chow ring also for non-complete varieties.

This is demonstrated by the following example. Note that the computation
for the non-complete variety leads to a Chow ring which is identical to
the Chow ring of a certain matroid. This observation can be anticipated
by e.g. the results in [FY04](@cite).

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> is_complete(p2)
true

julia> ngens(chow_ring(p2))
3

julia> v = normal_toric_variety(incidence_matrix([[1], [2], [3]]), [[1, 0], [0, 1], [-1, -1]])
Normal toric variety

julia> is_complete(v)
false

julia> set_coordinate_names(v, ["x_{1}", "x_{2}", "x_{3}"])

julia> chow_ring(v)
Quotient
  of multivariate polynomial ring in 3 variables x_{1}, x_{2}, x_{3}
    over rational field
  by ideal (x_{1} - x_{3}, x_{2} - x_{3}, x_{1}*x_{2}, x_{1}*x_{3}, x_{2}*x_{3})

julia> M = cycle_matroid(complete_graph(3))
Matroid of rank 2 on 3 elements

julia> chow_ring(M)
Quotient
  of multivariate polynomial ring in 3 variables x_{Edge(2, 1)}, x_{Edge(3, 1)}, x_{Edge(3, 2)}
    over rational field
  by ideal with 5 generators
```
"""
@attr MPolyQuoRing function chow_ring(v::NormalToricVarietyType)
    @req is_simplicial(v) "The combinatorial Chow ring is (currently) only supported for simplicial toric varieties"
    R, _ = polynomial_ring(coefficient_ring(v), coordinate_names(v); cached=false)
    linear_relations = ideal_of_linear_relations(R, v)
    stanley_reisner = stanley_reisner_ideal(R, v)
    return quo(R, linear_relations + stanley_reisner)[1]
end


@doc raw"""
    gens_of_rational_equivalence_classes(v::NormalToricVarietyType)

Return a list of generators of the Chow ring of a
complete, simplicial toric variety.

Recall that the cones of a complete, simplicial toric variety
can be seen as generators of the Chow ring (lemma 12.5.1 in
[CLS11](@cite)). This function first maps each cone to an
element of the Chow ring and then removes elements by taking
rational equivalence into account.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> gens_of_rational_equivalence_classes(p2)
6-element Vector{MPolyQuoRingElem{QQMPolyRingElem}}:
 x3^2
 x3^2
 x3^2
 x3
 x3
 x3
```
"""
@attr Vector{MPolyQuoRingElem{QQMPolyRingElem}} function gens_of_rational_equivalence_classes(v::NormalToricVarietyType)
  cr = chow_ring(v)
  R = base_ring(cr)
  cs = cones(v; trivial=false)
  return [simplify(cr(R([1], [Vector{Int}(cs[k,:])]))) for k in 1:nrows(cs)]
end


@doc raw"""
    map_gens_of_chow_ring_to_cox_ring(v::NormalToricVarietyType)

Return a dictionary which maps the generators of the chow
ring to monomials in the Cox ring. This dictionary involves
a choice, i.e. is not unique.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> p2 = projective_space(NormalToricVariety, 2);

julia> map_gens_of_chow_ring_to_cox_ring(p2)
Dict{QQMPolyRingElem, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 2 entries:
  x3^2 => x1*x3
  x3   => x3
```
"""
@attr Dict{QQMPolyRingElem, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} function map_gens_of_chow_ring_to_cox_ring(v::NormalToricVarietyType)
  cr = chow_ring(v)
  R = base_ring(cr)
  co = cox_ring(v)
  cs = cones(v; trivial=false)
  mapping = Dict{QQMPolyRingElem, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}()
  for k in 1:nrows(cs)
    p1 = simplify(cr(R([1], [Vector{Int}(cs[k,:])]))).f
    p2 = co([1], [Vector{Int}(cs[k,:])])
    coeff = [c for c in AbstractAlgebra.coefficients(p1)][1]
    if coeff != 1
      p1 = 1//coeff * p1
      p2 = 1//coeff * p2
    end
    mapping[p1] = p2
  end
  return mapping
end
