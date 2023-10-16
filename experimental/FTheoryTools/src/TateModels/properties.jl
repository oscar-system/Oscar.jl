#####################################################
# 1. Properties for general global Tate models
#####################################################

@doc raw"""
    base_fully_specified(t::GlobalTateModel)

Return `true` if the Tate model has a concrete base space and `false` otherwise.

```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> base_fully_specified(t)
false
```
"""
base_fully_specified(t::GlobalTateModel) = get_attribute(t, :base_fully_specified)