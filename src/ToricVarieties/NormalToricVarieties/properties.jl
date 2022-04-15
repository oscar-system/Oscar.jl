######################
# 4: Properties
######################

@doc Markdown.doc"""
    isnormal(v::AbstractNormalToricVariety)

Checks if the normal toric variety `v` is normal. (This function is somewhat tautological at this point.)

# Examples
```jldoctest
julia> isnormal(projective_space(NormalToricVariety, 2))
true
```
"""
function isnormal(v::AbstractNormalToricVariety)
    return true
end
export isnormal


@doc Markdown.doc"""
    isaffine(v::AbstractNormalToricVariety)

Checks if the normal toric variety `v` is affine.

# Examples
```jldoctest
julia> isaffine(projective_space(NormalToricVariety, 2))
false
```
"""
@attr Bool isaffine(v::AbstractNormalToricVariety) = pm_object(v).AFFINE
export isaffine


@doc Markdown.doc"""
    isprojective(v::AbstractNormalToricVariety)

Checks if the normal toric variety `v` is projective, i.e. if the fan of `v` is the the normal fan of a polytope.

# Examples
```jldoctest
julia> isprojective(projective_space(NormalToricVariety, 2))
true
```
"""
@attr Bool isprojective(v::AbstractNormalToricVariety) = pm_object(v).PROJECTIVE
export isprojective


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
    if issmooth(v) == false
        return false
    end
    if isprojective(v) == false
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
    issmooth(v::AbstractNormalToricVariety)

Checks if the normal toric variety `v` is smooth.

# Examples
```jldoctest
julia> issmooth(projective_space(NormalToricVariety, 2))
true
```
"""
@attr Bool issmooth(v::AbstractNormalToricVariety) = pm_object(v).SMOOTH
export issmooth


@doc Markdown.doc"""
    iscomplete(v::AbstractNormalToricVariety)

Checks if the normal toric variety `v` is complete.

# Examples
```jldoctest
julia> iscomplete(projective_space(NormalToricVariety, 2))
true
```
"""
@attr Bool iscomplete(v::AbstractNormalToricVariety) = pm_object(v).COMPLETE
export iscomplete


@doc Markdown.doc"""
    hastorusfactor(v::AbstractNormalToricVariety)

Checks if the normal toric variety `v` has a torus factor.

# Examples
```jldoctest
julia> hastorusfactor(projective_space(NormalToricVariety, 2))
false
```
"""
@attr Bool hastorusfactor(v::AbstractNormalToricVariety) = Polymake.common.rank(rays(v)) < ambient_dim(v)
export hastorusfactor


@doc Markdown.doc"""
    isorbifold(v::AbstractNormalToricVariety)

Checks if the normal toric variety `v` is an orbifold.

# Examples
```jldoctest
julia> isorbifold(projective_space(NormalToricVariety, 2))
true
```
"""
@attr Bool isorbifold(v::AbstractNormalToricVariety) = pm_object(v).SIMPLICIAL
export isorbifold


@doc Markdown.doc"""
    issimplicial(v::AbstractNormalToricVariety)

Checks if the normal toric variety `v` is simplicial. Hence, this function works just as `isorbifold`. It is implemented for user convenience.

# Examples
```jldoctest
julia> issimplicial(projective_space(NormalToricVariety, 2))
true
```
"""
issimplicial(v::AbstractNormalToricVariety) = isorbifold(v)
export issimplicial


@doc Markdown.doc"""
    isgorenstein(v::AbstractNormalToricVariety)

Checks if the normal toric variety `v` is Gorenstein.

# Examples
```jldoctest
julia> isgorenstein(projective_space(NormalToricVariety, 2))
true
```
"""
@attr Bool isgorenstein(v::AbstractNormalToricVariety) = pm_object(v).GORENSTEIN
export isgorenstein


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
    isfano(v::AbstractNormalToricVariety)

Checks if the normal toric variety `v` is fano.

# Examples
```jldoctest
julia> isfano(projective_space(NormalToricVariety, 2))
true
```
"""
@attr Bool isfano(v::AbstractNormalToricVariety) = pm_object(v).FANO
export isfano
