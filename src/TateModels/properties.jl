@doc Markdown.doc"""
    base_fully_specified(t::GlobalTateModel)

Return `true` is the Tate model has a concrete base space and `false` otherwise.

```jldoctest
julia> using Oscar

julia> test_space = hirzebruch_surface(2) * projective_space(NormalToricVariety,1)
A normal toric variety

julia> test_space1 = blowup_on_ith_minimal_torus_orbit(test_space,1,"e1")
A normal toric variety

julia> test_space2 = blowup_on_ith_minimal_torus_orbit(test_space1,1,"e2")
A normal toric variety

julia> base = blowup_on_ith_minimal_torus_orbit(test_space2,1,"e3")
A normal toric variety

julia> t = GlobalTateModel(base)
A global Tate model over a concrete base

julia> base_fully_specified(t)
true
```
"""
base_fully_specified(t::GlobalTateModel) = get_attribute(t, :base_fully_specified)
export base_fully_specified
