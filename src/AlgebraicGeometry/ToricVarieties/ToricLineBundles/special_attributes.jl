@doc raw"""
    vanishing_sets(variety::NormalToricVarietyType)

Compute the vanishing sets of an abstract toric variety `v` by use of the cohomCalg algorithm.
"""
@attr Vector{ToricVanishingSet} function vanishing_sets(variety::NormalToricVarietyType)
  @req (is_simplicial(variety) && is_projective(variety)) "the currently implemented cohomCalg algorithm only applies to toric varieties that are simplicial and projective"
  denominator_contributions = contributing_denominators(variety)
  vs = ToricVanishingSet[]
  for i in 1:length(denominator_contributions)
    list_of_polyhedra = Polyhedron{QQFieldElem}[
      turn_denominator_into_polyhedron(variety, m) for m in denominator_contributions[i]
    ]
    push!(vs, ToricVanishingSet(variety, list_of_polyhedra, [i - 1]))
  end
  return vs
end

@doc raw"""
    immaculate_line_bundles(variety::NormalToricVarietyType)

Compute all immaculate line bundles as a toric vanishing set by
intersecting the vanishing sets for all cohomology indices.

# Examples
```jldoctest
julia> dP1 = del_pezzo_surface(NormalToricVariety, 1)
Normal toric variety

julia> ilb = immaculate_line_bundles(dP1)
Toric vanishing set for cohomology indices [0, 1, 2]

julia> polyhedra(ilb)
4-element Vector{Polyhedron{QQFieldElem}}:
 Polyhedron in ambient dimension 2
 Polyhedron in ambient dimension 2
 Polyhedron in ambient dimension 2
 Polyhedron in ambient dimension 2

julia> print_constraints(polyhedra(ilb)[1])
-x_1 <= 0
-x_1 + x_2 <= 0

julia> print_constraints(polyhedra(ilb)[2])
-x_1 + x_2 <= 0
x_2 <= -2

julia> print_constraints(polyhedra(ilb)[3])
-x_2 <= -1
x_1 - x_2 <= -2

julia> print_constraints(polyhedra(ilb)[4])
x_1 - x_2 <= -2
x_1 <= -3
```
"""
@attr Any function immaculate_line_bundles(variety::NormalToricVarietyType)
  @req (is_simplicial(variety) && is_projective(variety)) "the currently implemented cohomCalg algorithm only applies to toric varieties that are simplicial and projective"
  denominator_contributions = reduce(vcat, contributing_denominators(variety))
  list_of_polyhedra = Polyhedron{QQFieldElem}[
    turn_denominator_into_polyhedron(variety, m) for m in denominator_contributions
  ]
  return ToricVanishingSet(variety, list_of_polyhedra, collect(0:dim(variety)))
end
