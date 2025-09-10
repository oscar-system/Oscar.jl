############################
# Attributes
############################

@doc raw"""
    toric_variety(tvs::ToricVanishingSet)

Return the toric variety of the vanishing set `tvs`.

# Examples
```jldoctest
julia> dP1 = del_pezzo_surface(NormalToricVariety, 1)
Normal toric variety

julia> vs = vanishing_sets(dP1)
3-element Vector{ToricVanishingSet}:
 Toric vanishing set for cohomology indices [0]
 Toric vanishing set for cohomology indices [1]
 Toric vanishing set for cohomology indices [2]

julia> toric_variety(vs[3])
Normal, simplicial, projective, 2-dimensional toric variety without torusfactor
```
"""
toric_variety(tvs::ToricVanishingSet) = tvs.toric_variety


@doc raw"""
    polyhedra(tvs::ToricVanishingSet)

Return the vector of the polyhedra whose complement defines the vanishing set `tvs`.

# Examples
```jldoctest
julia> dP1 = del_pezzo_surface(NormalToricVariety, 1)
Normal toric variety

julia> vs = vanishing_sets(dP1)
3-element Vector{ToricVanishingSet}:
 Toric vanishing set for cohomology indices [0]
 Toric vanishing set for cohomology indices [1]
 Toric vanishing set for cohomology indices [2]

julia> polyhedra(vs[3])
1-element Vector{Polyhedron{QQFieldElem}}:
 Polyhedron in ambient dimension 2
```
"""
polyhedra(tvs::ToricVanishingSet) = tvs.ps


@doc raw"""
    cohomology_indices(tvs::ToricVanishingSet)

Return the cohomology indices of the toric vanishing set `tvs`.

# Examples
```jldoctest
julia> dP1 = del_pezzo_surface(NormalToricVariety, 1)
Normal toric variety

julia> vs = vanishing_sets(dP1)
3-element Vector{ToricVanishingSet}:
 Toric vanishing set for cohomology indices [0]
 Toric vanishing set for cohomology indices [1]
 Toric vanishing set for cohomology indices [2]

julia> cohomology_indices(vs[3])
1-element Vector{Int64}:
 2
```
"""
cohomology_indices(tvs::ToricVanishingSet) = tvs.cis

