############################
# Attributes
############################

@doc raw"""
    toric_variety(c::ClosedSubvarietyOfToricVariety)

When constructing a closed subvariety, a toric variety
must be provided in which the closed subvariety is contained.
This method returns this initially provided toric supervariety.

Note however that perse, a closed subvariety can be contained
in different non-isomorphic toric varieties.

# Examples
```jldoctest
julia> f2 = hirzebruch_surface(NormalToricVariety, 2);

julia> (t1, x1, t2, x2) = gens(cox_ring(f2));

julia> c = closed_subvariety_of_toric_variety(f2, [t1])
Closed subvariety of a normal toric variety

julia> toric_variety(c) == f2
true
```
"""
@attr NormalToricVarietyType toric_variety(c::ClosedSubvarietyOfToricVariety) =
  c.toric_variety

@doc raw"""
    defining_ideal(c::ClosedSubvarietyOfToricVariety)

When constructing a closed subvariety, an ideal in the
Cox ring of a normal toric variety must be provided.
This method returns this initially provided ideal.

# Examples
```jldoctest
julia> f2 = hirzebruch_surface(NormalToricVariety, 2);

julia> (t1, x1, t2, x2) = gens(cox_ring(f2));

julia> c = closed_subvariety_of_toric_variety(f2, [t1])
Closed subvariety of a normal toric variety

julia> defining_ideal(c) == ideal([t1])
true
```
"""
@attr MPolyIdeal defining_ideal(c::ClosedSubvarietyOfToricVariety) = c.defining_ideal

@doc raw"""
    radical(c::ClosedSubvarietyOfToricVariety)

When constructing a closed subvariety, an ideal in the
Cox ring of a normal toric variety must be provided.
This method returns the radical of this initially provided
ideal.

# Examples
```jldoctest
julia> f2 = hirzebruch_surface(NormalToricVariety, 2);

julia> (t1, x1, t2, x2) = gens(cox_ring(f2));

julia> c = closed_subvariety_of_toric_variety(f2, [t1])
Closed subvariety of a normal toric variety

julia> radical(c) == ideal([t1])
true
```
"""
@attr MPolyIdeal function radical(c::ClosedSubvarietyOfToricVariety)
  return radical(defining_ideal(c))
end
