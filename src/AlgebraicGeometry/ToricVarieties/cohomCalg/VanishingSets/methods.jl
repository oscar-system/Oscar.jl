@doc raw"""
    contains(tvs::ToricVanishingSet, l::ToricLineBundle)

Check if the toric line bundle `l` is contained in the toric vanishing set `tvs`.

# Examples
```jldoctest
julia> dP1 = del_pezzo_surface(NormalToricVariety, 1)
Normal toric variety

julia> l = toric_line_bundle(dP1, [3, 2])
Toric line bundle on a normal toric variety

julia> all_cohomologies(l)
3-element Vector{ZZRingElem}:
 7
 0
 0

julia> vs = vanishing_sets(dP1)
3-element Vector{ToricVanishingSet}:
 Toric vanishing set for cohomology indices [0]
 Toric vanishing set for cohomology indices [1]
 Toric vanishing set for cohomology indices [2]

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
  class = divisor_class(toric_divisor_class(l)).coeff
  class = [class[1, i] for i in 1:ncols(class)]
  for p in polyhedra(tvs)
    if class in p
      return false
    end
  end
  return true
end
