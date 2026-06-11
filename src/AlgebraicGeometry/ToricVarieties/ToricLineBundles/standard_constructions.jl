########################
# 1: Special attributes of toric varieties
########################

@doc raw"""
    structure_sheaf(v::NormalToricVarietyType)

Construct the structure sheaf of a normal toric variety.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> structure_sheaf(v)
Toric line bundle on a normal toric variety
```
"""
@attr ToricLineBundle structure_sheaf(v::NormalToricVarietyType) = toric_line_bundle(
  v, zero(picard_group_with_map(v)[1])
)

@doc raw"""
    anticanonical_bundle(v::NormalToricVarietyType)

Construct the anticanonical bundle of a normal toric variety.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> anticanonical_bundle(v)
Toric line bundle on a normal toric variety
```
"""
@attr ToricLineBundle anticanonical_bundle(v::NormalToricVarietyType) = toric_line_bundle(
  v, anticanonical_divisor(v)
)

@doc raw"""
    canonical_bundle(v::NormalToricVarietyType)

Construct the canonical bundle of a normal toric variety.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> canonical_bundle(v)
Toric line bundle on a normal toric variety
```
"""
@attr ToricLineBundle canonical_bundle(v::NormalToricVarietyType) = inv(
  anticanonical_bundle(v)
)

@doc raw"""
    trivial_line_bundle(v::NormalToricVarietyType)

Construct the trivial line bundle on a normal toric variety.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> l = trivial_line_bundle(v)
Toric line bundle on a normal toric variety

julia> is_trivial(l)
true
```
"""
@attr ToricLineBundle trivial_line_bundle(v::NormalToricVarietyType) = toric_line_bundle(
  v, trivial_divisor(v)
)
