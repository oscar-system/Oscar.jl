@doc Markdown.doc"""
    underlying_scheme(X::ToricCoveredScheme)

For a toric covered scheme ``X``, this returns
the underlying scheme. In other words, by applying
this method, you obtain a scheme that has forgotten
its toric origin.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> toric_scheme = ToricCoveredScheme(P2)
Scheme of a toric variety with fan spanned by RayVector{fmpq}[[1, 0], [0, 1], [-1, -1]]

julia> underlying_scheme(toric_scheme)
covered scheme with 3 affine patches in its default covering
```
"""
underlying_scheme(X::ToricCoveredScheme) = X.X
export underlying_scheme


@doc Markdown.doc"""
    normal_toric_variety(X::ToricCoveredScheme)

For a toric covered scheme ``X``, this returns
the underlying normal toric variety.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> toric_scheme = ToricCoveredScheme(P2)
Scheme of a toric variety with fan spanned by RayVector{fmpq}[[1, 0], [0, 1], [-1, -1]]

julia> normal_toric_variety(toric_scheme)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor
```
"""
normal_toric_variety(X::ToricCoveredScheme) = X.ntv
export normal_toric_variety


@doc Markdown.doc"""
    fan(X::ToricCoveredScheme)

For a toric covered scheme ``X``, this returns
the fan of the underlying toric variety.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> toric_scheme = ToricCoveredScheme(P2)
Scheme of a toric variety with fan spanned by RayVector{fmpq}[[1, 0], [0, 1], [-1, -1]]

julia> fan(toric_scheme)
A polyhedral fan in ambient dimension 2
```
"""
fan(X::ToricCoveredScheme) = fan(normal_toric_variety(X))
export fan
