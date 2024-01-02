#####################################################
# 1. Properties for general global Tate models
#####################################################

@doc raw"""
    is_base_space_fully_specified(h::HypersurfaceModel)

Return `true` is the hypersurface model has a concrete base space and `false` otherwise.

```jldoctest
julia> h = hypersurface_model_over_projective_space(2)
Hypersurface model over a concrete base

julia> is_base_space_fully_specified(h)
true
```
"""
is_base_space_fully_specified(h::HypersurfaceModel) = !(typeof(h.base_space) <: FamilyOfSpaces)


@doc raw"""
    is_partially_resolved(h::HypersurfaceModel)

Return `true` if resolution techniques were applies to the hypersurface model,
thereby potentially resolving its singularities. Otherwise, return `false`.

```jldoctest
julia> h = hypersurface_model_over_projective_space(2)
Hypersurface model over a concrete base

julia> is_partially_resolved(h)
false
```
"""
is_partially_resolved(h::HypersurfaceModel) = get_attribute(h, :partially_resolved)
