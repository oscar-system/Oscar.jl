is_trivial(l::ToricLineBundle) = is_principal(toric_divisor(l))


@doc raw"""
    is_basepoint_free(l::ToricLineBundle)

Return `true` if the toric line bundle `l` is basepoint free and `false` otherwise.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal toric variety

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
Normal toric variety

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
Normal toric variety

julia> is_very_ample(toric_line_bundle(F4, [1,0]))
false
```
"""
is_very_ample(l::ToricLineBundle) = is_very_ample(toric_divisor(l))


@doc raw"""
    is_immaculate(l::ToricLineBundle)

Return `true` if all sheaf cohomologies of `l` are trivial and `false` otherwise.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal toric variety

julia> l = toric_line_bundle(F4, [1,0])
Toric line bundle on a normal toric variety

julia> is_immaculate(toric_line_bundle(F4, [1,0]))
false

julia> all_cohomologies(l)
3-element Vector{ZZRingElem}:
 2
 0
 0
```
"""
@attr Bool is_immaculate(l::ToricLineBundle) = all(is_zero, all_cohomologies(l))
