@doc Markdown.doc"""
    is_smooth(X::ToricCoveredScheme)

For a toric covered scheme ``X``, this method
return `true` if ``X`` is smooth and `false` otherwise.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> toric_scheme = ToricCoveredScheme(P2)
Scheme of a toric variety with fan spanned by RayVector{fmpq}[[1, 0], [0, 1], [-1, -1]]

julia> is_smooth(toric_scheme)
true
```
"""
is_smooth(X::ToricCoveredScheme) = is_smooth(normal_toric_variety(X))
export is_smooth
