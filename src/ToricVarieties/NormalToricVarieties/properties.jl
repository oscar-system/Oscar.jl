######################
# 4: Properties
######################

@doc Markdown.doc"""
    is_normal(v::AbstractNormalToricVariety)

Checks if the normal toric variety `v` is normal. (This function is somewhat tautological at this point.)

# Examples
```jldoctest
julia> is_normal(projective_space(NormalToricVariety, 2))
true
```
"""
function is_normal(v::AbstractNormalToricVariety)
    return true
end
export is_normal


@doc Markdown.doc"""
    is_affine(v::AbstractNormalToricVariety)

Checks if the normal toric variety `v` is affine.

# Examples
```jldoctest
julia> is_affine(projective_space(NormalToricVariety, 2))
false
```
"""
@attr Bool is_affine(v::AbstractNormalToricVariety) = pm_object(v).AFFINE
export is_affine


@doc Markdown.doc"""
    is_projective(v::AbstractNormalToricVariety)

Checks if the normal toric variety `v` is projective, i.e. if the fan of `v` is the the normal fan of a polytope.

# Examples
```jldoctest
julia> is_projective(projective_space(NormalToricVariety, 2))
true
```
"""
@attr Bool is_projective(v::AbstractNormalToricVariety) = pm_object(v).PROJECTIVE
export is_projective


@doc Markdown.doc"""
    is_projective_space(v::AbstractNormalToricVariety)

Decides if the normal toric varieties `v` is a projective space.

# Examples
```jldoctest
julia> F5 = hirzebruch_surface(5)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> is_projective_space(F5)
false

julia> is_projective_space(projective_space(NormalToricVariety, 2))
true
```
"""
@attr Bool function is_projective_space(v::AbstractNormalToricVariety)
    if is_smooth(v) == false
        return false
    end
    if is_projective(v) == false
        return false
    end
    if rank(class_group(v)) > 1
        return false
    end
    w = [[Int(x) for x in transpose(g.coeff)] for g in gens(class_group(v))]
    for g in gens(class_group(v))
        g = [Int(x) for x in g.coeff if !iszero(x)]
        if length(g) > 1
            return false
        end
        if g[1] != 1
            return false
        end
    end
    return irrelevant_ideal(v) == ideal(gens(cox_ring(v)))
end
export is_projective_space


@doc Markdown.doc"""
    is_smooth(v::AbstractNormalToricVariety)

Checks if the normal toric variety `v` is smooth.

# Examples
```jldoctest
julia> is_smooth(projective_space(NormalToricVariety, 2))
true
```
"""
@attr Bool is_smooth(v::AbstractNormalToricVariety) = pm_object(v).SMOOTH
export is_smooth


@doc Markdown.doc"""
    is_complete(v::AbstractNormalToricVariety)

Checks if the normal toric variety `v` is complete.

# Examples
```jldoctest
julia> is_complete(projective_space(NormalToricVariety, 2))
true
```
"""
@attr Bool is_complete(v::AbstractNormalToricVariety) = pm_object(v).COMPLETE
export is_complete


@doc Markdown.doc"""
    has_torusfactor(v::AbstractNormalToricVariety)

Checks if the normal toric variety `v` has a torus factor.

# Examples
```jldoctest
julia> has_torusfactor(projective_space(NormalToricVariety, 2))
false
```
"""
@attr Bool has_torusfactor(v::AbstractNormalToricVariety) = Polymake.common.rank(rays(v)) < ambient_dim(v)
export has_torusfactor


@doc Markdown.doc"""
    is_orbifold(v::AbstractNormalToricVariety)

Checks if the normal toric variety `v` is an orbifold.

# Examples
```jldoctest
julia> is_orbifold(projective_space(NormalToricVariety, 2))
true
```
"""
@attr Bool is_orbifold(v::AbstractNormalToricVariety) = pm_object(v).SIMPLICIAL
export is_orbifold


@doc Markdown.doc"""
    is_simplicial(v::AbstractNormalToricVariety)

Checks if the normal toric variety `v` is simplicial. Hence, this function works just as `is_orbifold`. It is implemented for user convenience.

# Examples
```jldoctest
julia> is_simplicial(projective_space(NormalToricVariety, 2))
true
```
"""
is_simplicial(v::AbstractNormalToricVariety) = is_orbifold(v)
export is_simplicial


@doc Markdown.doc"""
    is_gorenstein(v::AbstractNormalToricVariety)

Checks if the normal toric variety `v` is Gorenstein.

# Examples
```jldoctest
julia> is_gorenstein(projective_space(NormalToricVariety, 2))
true
```
"""
@attr Bool is_gorenstein(v::AbstractNormalToricVariety) = pm_object(v).GORENSTEIN
export is_gorenstein


@doc Markdown.doc"""
    is_q_gorenstein(v::AbstractNormalToricVariety)

Checks if the normal toric variety `v` is Q-Gorenstein.

# Examples
```jldoctest
julia> is_q_gorenstein(projective_space(NormalToricVariety, 2))
true
```
"""
@attr Bool is_q_gorenstein(v::AbstractNormalToricVariety) = pm_object(v).Q_GORENSTEIN
export is_q_gorenstein


@doc Markdown.doc"""
    is_fano(v::AbstractNormalToricVariety)

Checks if the normal toric variety `v` is fano.

# Examples
```jldoctest
julia> is_fano(projective_space(NormalToricVariety, 2))
true
```
"""
@attr Bool is_fano(v::AbstractNormalToricVariety) = pm_object(v).FANO
export is_fano
