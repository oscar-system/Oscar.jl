@doc Markdown.doc"""
    base_fully_specified(t::GlobalTateModel)

Return `true` is the Tate model has a concrete base space and `false` otherwise.

```jldoctest
julia> t = su5_tate_model_over_arbitrary_3d_base()
Global Tate model over a not fully specified base

julia> base_fully_specified(t)
false
```
"""
base_fully_specified(t::GlobalTateModel) = get_attribute(t, :base_fully_specified)
