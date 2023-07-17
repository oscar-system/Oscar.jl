@doc raw"""
    base_fully_specified(w::WeierstrassModel)

Return `true` is the Weierstrass model has a concrete base space and `false` otherwise.

```jldoctest
julia> w = su5_weierstrass_model_over_arbitrary_3d_base()
Assuming that the first row of the given grading is the grading under Kbar

Weierstrass model over a not fully specified base

julia> base_fully_specified(w)
false
```
"""
base_fully_specified(w::WeierstrassModel) = get_attribute(w, :base_fully_specified)
