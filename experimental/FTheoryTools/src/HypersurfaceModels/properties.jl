#####################################################
# 1. Properties for general global Tate models
#####################################################

@doc raw"""
    base_fully_specified(h::HypersurfaceModel)

Return `true` is the hypersurface model has a concrete base space and `false` otherwise.

```jldoctest
julia> h = hypersurface_model_over_projective_space(2)
Hypersurface model over a concrete base

julia> base_fully_specified(h)
true
```
"""
base_fully_specified(h::HypersurfaceModel) = get_attribute(h, :base_fully_specified)
