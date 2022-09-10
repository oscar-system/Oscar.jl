######################
# Properties
######################

@doc Markdown.doc"""
    is_empty(c::ClosedSubvarietyOfToricVariety)

Checks if a closed subvariety of a toric variety is empty.
This check uses proposition 5.2.6 in [CLS11](@cite).

# Examples
```jldoctest
julia> f2 = hirzebruch_surface(2);

julia> (t1, x1, t2, x2) = gens(cox_ring(f2));

julia> c = ClosedSubvarietyOfToricVariety(f2, [t1])
A closed subvariety of a normal toric variety

julia> is_empty(c)
false
```
"""
@attr Bool function is_empty(c::ClosedSubvarietyOfToricVariety)
    B = irrelevant_ideal(toric_variety(c))
    return intersect(radical(c), B) == B
end
export is_empty
