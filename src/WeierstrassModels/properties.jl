@doc Markdown.doc"""
    base_fully_specified(t::GlobalWeierstrassModel)

Return `true` is the Weierstrass model has a concrete base space and `false` otherwise.

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

julia> w = GlobalWeierstrassModel(base)
A global Weierstrass model over a concrete base

julia> base_fully_specified(w)
true
```
"""
base_fully_specified(w::GlobalWeierstrassModel) = get_attribute(w, :base_fully_specified)
export base_fully_specified
