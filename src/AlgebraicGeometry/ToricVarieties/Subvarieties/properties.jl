######################
# Properties
######################

@doc raw"""
    is_empty(c::ClosedSubvarietyOfToricVariety)

Check if a closed subvariety of a toric variety is empty.
This check uses proposition 5.2.6 in [CLS11](@cite).

# Examples
```jldoctest
julia> f2 = hirzebruch_surface(NormalToricVariety, 2);

julia> (t1, x1, t2, x2) = gens(cox_ring(f2));

julia> c = closed_subvariety_of_toric_variety(f2, [t1])
Closed subvariety of a normal toric variety

julia> is_empty(c)
false

julia> c2 = closed_subvariety_of_toric_variety(f2, [x1,x2])
Closed subvariety of a normal toric variety

julia> is_empty(c2)
true
```
"""
@attr Bool function is_empty(c::ClosedSubvarietyOfToricVariety)
  B = irrelevant_ideal(toric_variety(c))
  return (is_one(defining_ideal(c)) || is_one(saturation(defining_ideal(c), B)))
end
