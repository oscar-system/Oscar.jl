@doc Markdown.doc"""
    base_fully_specified(t::GlobalWeierstrassModel)

Return `true` is the Weierstrass model has a concrete base space and `false` otherwise.

```jldoctest
julia> auxiliary_base_ring, (f, g, x) = QQ["f", "g", "x"];

julia> w = global_weierstrass_model(f, g, auxiliary_base_ring, 3)
Global Weierstrass model over a not fully specified base

julia> base_fully_specified(w)
false
```
"""
base_fully_specified(w::GlobalWeierstrassModel) = get_attribute(w, :base_fully_specified)
