############################
# Attributes
############################

@doc Markdown.doc"""
    toric_variety(tvs::ToricVanishingSet)

Return the toric variety of the vanishing set `tvs`.

# Examples
```jldoctest
julia> dP1 = del_pezzo(1)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> vs = vanishing_sets(dP1)
3-element Vector{ToricVanishingSet}:
 A toric vanishing set for cohomology index 0
 A toric vanishing set for cohomology index 1
 A toric vanishing set for cohomology index 2

julia> toric_variety(vs[3])
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety over QQ without torusfactor
```
"""
function toric_variety(tvs::ToricVanishingSet)
    return tvs.toric_variety::AbstractNormalToricVariety
end
export toric_variety


@doc Markdown.doc"""
    polyhedra(tvs::ToricVanishingSet)

Return the vector of the polyhedra whose complement defines the vanishing set `tvs`.

# Examples
```jldoctest
julia> dP1 = del_pezzo(1)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> vs = vanishing_sets(dP1)
3-element Vector{ToricVanishingSet}:
 A toric vanishing set for cohomology index 0
 A toric vanishing set for cohomology index 1
 A toric vanishing set for cohomology index 2

julia> polyhedra(vs[3])
1-element Vector{Polyhedra{fmpq}}:
 A polyhedron in ambient dimension 2

"""
function polyhedra(tvs::ToricVanishingSet)
    return tvs.ps::Vector{Polyhedron{fmpq}}
end
export polyhedra


@doc Markdown.doc"""
    cohomology_index(tvs::ToricVanishingSet)

Return the cohomology index of the toric vanishing set `tvs`.

# Examples
```jldoctest
julia> dP1 = del_pezzo(1)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> vs = vanishing_sets(dP1)
3-element Vector{ToricVanishingSet}:
 A toric vanishing set for cohomology index 0
 A toric vanishing set for cohomology index 1
 A toric vanishing set for cohomology index 2

julia> cohomology_index(vs[3])
2
```
"""
function cohomology_index(tvs::ToricVanishingSet)
    return tvs.i::Int
end
export cohomology_index
