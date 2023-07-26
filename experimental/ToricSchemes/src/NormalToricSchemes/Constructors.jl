###############################
### 1: General constructor
###############################

@doc raw"""
    toric_spec(antv::AffineNormalToricVariety)

Constructs the affine toric scheme (i.e. Spec of a ring) associated to an affine toric variety.

# Examples
```jldoctest
julia> C = positive_hull([1 0; 0 1])
Polyhedral cone in ambient dimension 2

julia> antv = affine_normal_toric_variety(C)
Normal, affine toric variety

julia> toric_spec(antv)
Spec of an affine toric variety
```
"""
toric_spec(antv::AffineNormalToricVariety) = ToricSpec(antv)


@doc raw"""
    toric_covered_scheme(antv::AffineNormalToricVariety)

Constructs the toric covered scheme associated to an affine toric variety.

# Examples
```jldoctest
julia> C = positive_hull([1 0; 0 1])
Polyhedral cone in ambient dimension 2

julia> antv = affine_normal_toric_variety(C)
Normal, affine toric variety

julia> toric_covered_scheme(antv)
Scheme of a toric variety
```
"""
toric_covered_scheme(antv::AffineNormalToricVariety) = ToricCoveredScheme(normal_toric_variety(polyhedral_fan(antv)))


@doc raw"""
    toric_covered_scheme(ntv::NormalToricVariety)

Constructs the toric covered scheme associated to a normal toric variety.

# Examples
```jldoctest
julia> ntv = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> toric_covered_scheme(ntv)
Scheme of a toric variety
```
"""
toric_covered_scheme(ntv::NormalToricVariety) = ToricCoveredScheme(ntv)


###############################
### 2: Special constructors
###############################

@doc raw"""
    affine_space(::Type{ToricCoveredScheme}, d::Int; set_attributes::Bool = true)

Constructs the affine space of dimension `d` as toric (covered) scheme.

# Examples
```jldoctest
julia> affine_space(ToricCoveredScheme, 2)
Scheme of a toric variety
```
"""
affine_space(::Type{ToricCoveredScheme}, d::Int; set_attributes::Bool = true) = ToricCoveredScheme(affine_space(NormalToricVariety, d; set_attributes = set_attributes))


@doc raw"""
    projective_space(::Type{ToricCoveredScheme}, d::Int; set_attributes::Bool = true)

Construct the projective space of dimension `d` as toric (covered) scheme.

# Examples
```jldoctest
julia> projective_space(ToricCoveredScheme, 2)
Scheme of a toric variety
```
"""
projective_space(::Type{ToricCoveredScheme}, d::Int; set_attributes::Bool = true) = ToricCoveredScheme(projective_space(NormalToricVariety, d; set_attributes = set_attributes))


@doc raw"""
    weighted_projective_space(::Type{ToricCoveredScheme}, w::Vector{T}; set_attributes::Bool = true) where {T <: IntegerUnion}

Construct the weighted projective space corresponding to the weights `w` as toric (covered) scheme.

# Examples
```jldoctest
julia> weighted_projective_space(ToricCoveredScheme, [2,3,1])
Scheme of a toric variety
```
"""
weighted_projective_space(::Type{ToricCoveredScheme}, w::Vector{T}; set_attributes::Bool = true) where {T <: IntegerUnion} = ToricCoveredScheme(weighted_projective_space(NormalToricVariety, w; set_attributes = set_attributes))


@doc raw"""
    hirzebruch_surface(::Type{ToricCoveredScheme}, r::Int; set_attributes::Bool = true)

Constructs the r-th Hirzebruch surface as toric (covered) scheme.

# Examples
```jldoctest
julia> hirzebruch_surface(ToricCoveredScheme, 5)
Scheme of a toric variety
```
"""
hirzebruch_surface(::Type{ToricCoveredScheme}, r::Int; set_attributes::Bool = true) = ToricCoveredScheme(hirzebruch_surface(NormalToricVariety, r; set_attributes = set_attributes))


@doc raw"""
    del_pezzo_surface(::Type{ToricCoveredScheme}, b::Int; set_attributes::Bool = true)

Constructs the del Pezzo surface with `b` blowups for `b` at most 3 as toric (covered) scheme.

# Examples
```jldoctest
julia> del_pezzo_surface(ToricCoveredScheme, 3)
Scheme of a toric variety
```
"""
del_pezzo_surface(::Type{ToricCoveredScheme}, b::Int; set_attributes::Bool = true) = ToricCoveredScheme(del_pezzo_surface(NormalToricVariety, b; set_attributes = set_attributes))
