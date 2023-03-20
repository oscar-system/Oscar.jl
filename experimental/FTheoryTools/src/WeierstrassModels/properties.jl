@doc Markdown.doc"""
    base_fully_specified(t::GlobalWeierstrassModel)

Return `true` is the Weierstrass model has a concrete base space and `false` otherwise.

```jldoctest
julia> using Oscar

julia> w = global_weierstrass_model(test_base())
A global Weierstrass model over a concrete base

julia> base_fully_specified(w)
true
```
"""
base_fully_specified(w::GlobalWeierstrassModel) = get_attribute(w, :base_fully_specified)
export base_fully_specified
