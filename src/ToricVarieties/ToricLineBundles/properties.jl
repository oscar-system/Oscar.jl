@doc Markdown.doc"""
    istrivial(l::ToricLineBundle)

Returns "true" if the toric line bundle `l` is trivial and "false" otherwise.

# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> istrivial(ToricLineBundle(H, [1,0]))
false

julia> istrivial(ToricLineBundle(H, [0,0]))
true
```
"""
istrivial(l::ToricLineBundle) = isprincipal(toric_divisor(l))
export istrivial


@doc Markdown.doc"""
    is_basepoint_free(l::ToricLineBundle)

Returns "true" if the toric line bundle `l` is basepoint free and "false" otherwise.

# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> is_basepoint_free(ToricLineBundle(H, [1,0]))
true
```
"""
is_basepoint_free(l::ToricLineBundle) = is_basepoint_free(toric_divisor(l))
export is_basepoint_free


@doc Markdown.doc"""
    isample(l::ToricLineBundle)

Returns "true" if the toric line bundle `l` is ample and "false" otherwise.

# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> isample(ToricLineBundle(H, [1,0]))
false
```
"""
isample(l::ToricLineBundle) = isample(toric_divisor(l))
export isample


@doc Markdown.doc"""
    is_very_ample(l::ToricLineBundle)

Returns "true" if the toric line bundle `l` is very ample and "false" otherwise.

# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> is_very_ample(ToricLineBundle(H, [1,0]))
false
```
"""
is_very_ample(l::ToricLineBundle) = is_very_ample(toric_divisor(l))
export is_very_ample
