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


@doc raw"""
    is_partially_resolved(w::WeierstrassModel)

Return `true` if resolution techniques were applies to the Weierstrass model,
thereby potentially resolving its singularities.

```jldoctest
julia> w = literature_model(arxiv_id = "1208.2695", equation = "B.19", completeness_check = false)
Assuming that the first row of the given grading is the grading under Kbar

Weierstrass model over a not fully specified base -- U(1) Weierstrass model based on arXiv paper 1208.2695 Eq. (B.19)

julia> is_partially_resolved(w)
false
```
"""
is_partially_resolved(w::WeierstrassModel) = get_attribute(w, :partially_resolved)