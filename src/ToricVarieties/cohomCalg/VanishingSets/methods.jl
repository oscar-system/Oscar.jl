@doc Markdown.doc"""
    contains(tvs::ToricVanishingSet, l::ToricLineBundle)

Checks if the toric line bundle `l` is contained in the toric vanishing set `tvs`.

# Examples
```jldoctest
julia> dP1 = del_pezzo_surface(1)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> l = ToricLineBundle(dP1, [3, 2])
A toric line bundle on a normal toric variety

julia> all_cohomologies(l)
3-element Vector{fmpz}:
 7
 0
 0

julia> vs = vanishing_sets(dP1)
3-element Vector{ToricVanishingSet}:
 A toric vanishing set for cohomology index 0
 A toric vanishing set for cohomology index 1
 A toric vanishing set for cohomology index 2

julia> contains(vs[1], l)
false

julia> contains(vs[2], l)
true

julia> contains(vs[3], l)
true
```
"""
function contains(tvs::ToricVanishingSet, l::ToricLineBundle)
    if toric_variety(l) !== toric_variety(tvs)
        return false
    end
    class = divisor_class(l).coeff
    class = [class[1, i] for i in 1:ncols(class)]
    for p in polyhedra(tvs)
        if contains(p, class)
            return false
        end
    end
    return true
end
export contains