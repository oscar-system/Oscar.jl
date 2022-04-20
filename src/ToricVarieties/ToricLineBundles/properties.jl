is_trivial(l::ToricLineBundle) = is_principal(toric_divisor(l))
export istrivial


@doc Markdown.doc"""
    is_basepoint_free(l::ToricLineBundle)

Return `true` if the toric line bundle `l` is basepoint free and `false` otherwise.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> is_basepoint_free(ToricLineBundle(F4, [0,1]))
true
```
"""
is_basepoint_free(l::ToricLineBundle) = is_basepoint_free(toric_divisor(l))
export is_basepoint_free


@doc Markdown.doc"""
    is_ample(l::ToricLineBundle)

Return `true` if the toric line bundle `l` is ample and `false` otherwise.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> is_ample(ToricLineBundle(F4, [1,0]))
false
```
"""
is_ample(l::ToricLineBundle) = is_ample(toric_divisor(l))
export is_ample


@doc Markdown.doc"""
    is_very_ample(l::ToricLineBundle)

Return `true` if the toric line bundle `l` is very ample and `false` otherwise.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> is_very_ample(ToricLineBundle(F4, [1,0]))
false
```
"""
is_very_ample(l::ToricLineBundle) = is_very_ample(toric_divisor(l))
export is_very_ample
