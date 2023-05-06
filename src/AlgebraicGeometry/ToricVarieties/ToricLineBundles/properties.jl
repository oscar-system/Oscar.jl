is_trivial(l::ToricLineBundle) = is_principal(toric_divisor(l))


@doc raw"""
    is_basepoint_free(l::ToricLineBundle)

Return `true` if the toric line bundle `l` is basepoint free and `false` otherwise.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> is_basepoint_free(toric_line_bundle(F4, [1, 0]))
true
```
"""
is_basepoint_free(l::ToricLineBundle) = is_basepoint_free(toric_divisor(l))


@doc raw"""
    is_ample(l::ToricLineBundle)

Return `true` if the toric line bundle `l` is ample and `false` otherwise.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> is_ample(toric_line_bundle(F4, [1,0]))
false
```
"""
is_ample(l::ToricLineBundle) = is_ample(toric_divisor(l))


@doc raw"""
    is_very_ample(l::ToricLineBundle)

Return `true` if the toric line bundle `l` is very ample and `false` otherwise.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> is_very_ample(toric_line_bundle(F4, [1,0]))
false
```
"""
is_very_ample(l::ToricLineBundle) = is_very_ample(toric_divisor(l))
