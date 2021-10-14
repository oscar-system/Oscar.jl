######################
# 4: Properties
######################

@doc Markdown.doc"""
    isnormal( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` is normal. (This function is somewhat tautological at this point.)

# Examples
```jldoctest
julia> isnormal(projective_space( 2 ))
true
```
"""
function isnormal( v::AbstractNormalToricVariety )
    return true
end
export isnormal


@doc Markdown.doc"""
    isaffine( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` is affine.

# Examples
```jldoctest
julia> isaffine( projective_space( 2 ) )
false
```
"""
function isaffine( v::AbstractNormalToricVariety )
    return pm_ntv(v).AFFINE::Bool
end
export isaffine


@doc Markdown.doc"""
    isprojective( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` is projective, i.e. if the fan of `v` is the the normal fan of a polytope.

# Examples
```jldoctest
julia> isprojective( projective_space( 2 ) )
true
```
"""
function isprojective( v::AbstractNormalToricVariety )
    return pm_ntv(v).PROJECTIVE::Bool
end
export isprojective


@doc Markdown.doc"""
    issmooth( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` is smooth.

# Examples
```jldoctest
julia> issmooth( projective_space( 2 ) )
true
```
"""
function issmooth( v::AbstractNormalToricVariety )
    return pm_ntv(v).SMOOTH::Bool
end
export issmooth


@doc Markdown.doc"""
    iscomplete( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` is complete.

# Examples
```jldoctest
julia> iscomplete( projective_space( 2 ) )
true
```
"""
function iscomplete( v::AbstractNormalToricVariety )
    return pm_ntv(v).COMPLETE::Bool
end
export iscomplete


@doc Markdown.doc"""
    hastorusfactor( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` has a torus factor.

# Examples
```jldoctest
julia> hastorusfactor( projective_space( 2 ) )
false
```
"""
function hastorusfactor( v::AbstractNormalToricVariety )
    return v.polymakeNTV.FAN_DIM < v.polymakeNTV.FAN_AMBIENT_DIM
end
export hastorusfactor


@doc Markdown.doc"""
    isorbifold( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` is an orbifold.

# Examples
```jldoctest
julia> isorbifold( projective_space( 2 ) )
true
```
"""
function isorbifold( v::AbstractNormalToricVariety )
    return pm_ntv(v).SIMPLICIAL::Bool
end
export isorbifold


@doc Markdown.doc"""
    issimplicial( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` is simplicial. Hence, this function works just as `isorbifold`. It is implemented for user convenience.

# Examples
```jldoctest
julia> issimplicial( projective_space( 2 ) )
true
```
"""
function issimplicial( v::AbstractNormalToricVariety )
    return pm_ntv(v).SIMPLICIAL::Bool
end
export issimplicial


@doc Markdown.doc"""
    isgorenstein( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` is Gorenstein.

# Examples
```jldoctest
julia> isgorenstein( projective_space( 2 ) )
true
```
"""
function isgorenstein( v::AbstractNormalToricVariety )
    return pm_ntv(v).GORENSTEIN::Bool
end
export isgorenstein


@doc Markdown.doc"""
    isq_gorenstein( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` is Q-Gorenstein.

# Examples
```jldoctest
julia> isq_gorenstein( projective_space( 2 ) )
true
```
"""
function isq_gorenstein( v::AbstractNormalToricVariety )
    return pm_ntv(v).Q_GORENSTEIN::Bool
end
export isq_gorenstein


@doc Markdown.doc"""
    isfano( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` is fano.

# Examples
```jldoctest
julia> isfano( projective_space( 2 ) )
true
```
"""
function isfano( v::AbstractNormalToricVariety )
    return pm_ntv(v).FANO::Bool
end
export isfano
