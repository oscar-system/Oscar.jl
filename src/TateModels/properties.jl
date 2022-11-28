@doc Markdown.doc"""
    base_fully_specified(t::GlobalTateModel)

Return `true` is the Tate model has a concrete base space and `false` otherwise.

```jldoctest
julia> using Oscar

julia> t = GlobalTateModel(TestBase())
A global Tate model over a concrete base

julia> base_fully_specified(t)
true
```
"""
base_fully_specified(t::GlobalTateModel) = get_attribute(t, :base_fully_specified)
export base_fully_specified
