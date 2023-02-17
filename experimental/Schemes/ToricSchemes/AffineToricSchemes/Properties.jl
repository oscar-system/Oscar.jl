@doc Markdown.doc"""
    is_smooth(X::ToricSpec)

Checks whether this toric scheme is smooth.
If so, the function returns `true` and otherwise
`false`.

# Examples
```jldoctest
julia> C = positive_hull([1 0; 0 1])
A polyhedral cone in ambient dimension 2

julia> antv = AffineNormalToricVariety(C)
A normal, affine toric variety

julia> affine_toric_scheme = ToricSpec(antv)
Spec of an affine toric variety with cone spanned by RayVector{fmpq}[[1, 0], [0, 1]]

julia> is_smooth(affine_toric_scheme)
true

julia> C2 = positive_hull([-1 1; 1 1])
A polyhedral cone in ambient dimension 2

julia> antv2 = AffineNormalToricVariety(C2)
A normal, affine toric variety

julia> affine_toric_scheme2 = ToricSpec(antv2)
Spec of an affine toric variety with cone spanned by RayVector{fmpq}[[-1, 1], [1, 1]]

julia> is_smooth(affine_toric_scheme2)
false
```
"""
is_smooth(X::ToricSpec) = is_smooth(affine_normal_toric_variety(X))
export is_smooth


@doc Markdown.doc"""
    is_simplicial(X::ToricSpec)

Checks whether this toric scheme is simplicial.
If so, the function returns `true` and otherwise
`false`.

# Examples
```jldoctest
julia> C = positive_hull([-1 1; 1 1])
A polyhedral cone in ambient dimension 2

julia> antv = AffineNormalToricVariety(C)
A normal, affine toric variety

julia> affine_toric_scheme = ToricSpec(antv)
Spec of an affine toric variety with cone spanned by RayVector{fmpq}[[-1, 1], [1, 1]]

julia> is_smooth(affine_toric_scheme)
false

julia> is_smooth(underlying_scheme(affine_toric_scheme))
false

julia> is_simplicial(affine_toric_scheme)
true
```
"""
is_simplicial(X::ToricSpec) = is_simplicial(affine_normal_toric_variety(X))
export is_simplicial
