######################
# 1: The Julia type for ToricVarieties
######################
abstract type AbstractNormalToricVariety end

struct NormalToricVariety <: AbstractNormalToricVariety
           polymakeNTV::Polymake.BigObject
end
export NormalToricVariety

struct AffineNormalToricVariety <: AbstractNormalToricVariety
           polymakeNTV::Polymake.BigObject
end
export AffineNormalToricVariety


function pm_ntv(v::AbstractNormalToricVariety)
    return v.polymakeNTV
end

######################
# 2: Generic constructors
######################


@doc Markdown.doc"""
    AffineNormalToricVariety(C::Cone)

Construct the affine normal toric variety $U_{C}$ corresponding to a polyhedral
cone `C`.

# Examples
Set `C` to be the positive orthant in two dimensions.
```jldoctest
julia> C = Oscar.positive_hull([1 0; 0 1])
A polyhedral cone in ambient dimension 2

julia> antv = AffineNormalToricVariety(C)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""
function AffineNormalToricVariety(C::Cone)
    pmc = Oscar.pm_cone(C)
    fan = Polymake.fan.check_fan_objects(pmc)
    pmntv = Polymake.fulton.NormalToricVariety(fan)
    return AffineNormalToricVariety(pmntv)
end


@doc Markdown.doc"""
    NormalToricVariety(C::Cone)
Construct the (affine) normal toric variety $X_{\Sigma}$ corresponding to a
polyhedral fan $\Sigma = C$ consisting only of the cone `C`.
# Examples
Set `C` to be the positive orthant in two dimensions.
```jldoctest
julia> C = Oscar.positive_hull([1 0; 0 1])
A polyhedral cone in ambient dimension 2
julia> ntv = NormalToricVariety(C)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""
function NormalToricVariety(C::Cone)
    return AffineNormalToricVariety(C)
end


@doc Markdown.doc"""
    NormalToricVariety(PF::PolyhedralFan)

Construct the normal toric variety $X_{PF}$ corresponding to a polyhedral fan `PF`.

# Examples
Take `PF` to be the normal fan of the square.
```jldoctest
julia> square = Oscar.cube(2)
A polyhedron in ambient dimension 2

julia> nf = Oscar.normal_fan(square)
A polyhedral fan in ambient dimension 2

julia> ntv = NormalToricVariety(nf)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""    
function NormalToricVariety(PF::PolyhedralFan)
    fan = Oscar.pm_fan(PF)
    pmntv = Polymake.fulton.NormalToricVariety(fan)
    if fan.N_MAXIMAL_CONES == 1
        return AffineNormalToricVariety( pmntv )
    end
    return NormalToricVariety(pmntv)
end


@doc Markdown.doc"""
    NormalToricVariety(P::Polyhedron)

Construct the normal toric variety $X_{\Sigma_P}$ corresponding to the normal
fan $\Sigma_P$ of the given polyhedron `P`.

Note that this only coincides with the projective variety associated to `P`, if
`P` is normal.

# Examples
Set `P` to be a square.
```jldoctest
julia> square = Oscar.cube(2)
A polyhedron in ambient dimension 2

julia> ntv = NormalToricVariety(square)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""    
function NormalToricVariety(P::Polyhedron)
    fan = normal_fan(P)
    return NormalToricVariety(fan)
end


@doc Markdown.doc"""
    NormalToricVariety( r::Matrix{Int}, c::Vector{Vector{Int}} )

Construct the normal toric variety whose fan has ray generators `r` and maximal cones `c`.

# Examples
```jldoctest
julia> NormalToricVariety( [-1 5; 0 1; 1 0; 0 -1], [[1,2],[2,3],[3,4],[4,1]] )
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""
function NormalToricVariety( rays::Matrix{Int}, cones::Vector{Vector{Int}} )
    Incidence = Oscar.IncidenceMatrix(cones)
    arr = Polymake.@convert_to Array{Set{Int}} Polymake.common.rows(Incidence.pm_incidencematrix)
    pmntv = Polymake.fulton.NormalToricVariety(
        RAYS = Oscar.matrix_for_polymake(rays),
        MAXIMAL_CONES = arr,
    )
    if length( cones ) == 1
        return AffineNormalToricVariety( pmntv )
    end    
    return NormalToricVariety( pmntv )
end

export NormalToricVariety



######################
# 3: Special constructors
######################

@doc Markdown.doc"""
    projective_space( d::Int )

Construct the projective space of dimension `d`.

# Examples
```jldoctest
julia> projective_space( 2 )
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""
function projective_space( d::Int )
    f = normal_fan(Oscar.simplex(d))
    pm_ntv = Polymake.fulton.NormalToricVariety(Oscar.pm_fan(f))
    return NormalToricVariety(pm_ntv)
end
export projective_space


@doc Markdown.doc"""
    hirzebruch_surface( r::Int )

Constructs the r-th Hirzebruch surface.

# Examples
```jldoctest
julia> hirzebruch_surface( 5 )
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""
function hirzebruch_surface( r::Int )
    Rays = [ 1 0; 0 1; -1 r; 0 -1]
    Cones = [[1,2],[2,3],[3,4],[4,1]]
    return NormalToricVariety( Rays, Cones )
end
export hirzebruch_surface


@doc Markdown.doc"""
    delPezzo( b::Int )

Constructs the delPezzo surface with b blowups for b at most 3.

# Examples
```jldoctest
julia> del_pezzo( 3 )
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""
function del_pezzo( b::Int )
    if b < 0
        throw(ArgumentError("Number of blowups for construction of delPezzo surfaces must be non-negative."))
        return 0
    end
    if b == 0 
        return projective_space( 2 )
    end
    if b == 1
        Rays = [ 1 0; 0 1; -1 0; -1 -1 ]
        Cones = [ [1,2],[2,3],[3,4],[4,1] ]
        return NormalToricVariety( Rays, Cones )
    end
    if b == 2
        Rays = [ 1 0; 0 1; -1 0; -1 -1; 0 -1 ]
        Cones = [ [1,2],[2,3],[3,4],[4,5],[5,1] ]
        return NormalToricVariety( Rays, Cones )
    end
    if b == 3
        Rays = [ 1 0; 1 1; 0 1; -1 0; -1 -1; 0 -1 ]
        Cones = [ [1,2],[2,3],[3,4],[4,5],[5,6],[6,1] ]
        return NormalToricVariety( Rays, Cones )
    end
    if b > 3
        throw(ArgumentError("delPezzo surfaces with more than 3 blowups are realized as subvarieties of toric ambient spaces. This is currently not supported."))
        return 0
    end
end
export del_pezzo


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
    has_torusfactor( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` has a torus factor.

# Examples
```jldoctest
julia> has_torusfactor( projective_space( 2 ) )
false
```
"""
function has_torusfactor( v::AbstractNormalToricVariety )
    return ( v.polymakeNTV.FAN_DIM < v.polymakeNTV.FAN_AMBIENT_DIM )::Bool
end
export has_torusfactor


@doc Markdown.doc"""
    is_orbifold( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` is an orbifold.

# Examples
```jldoctest
julia> is_orbifold( projective_space( 2 ) )
true
```
"""
function is_orbifold( v::AbstractNormalToricVariety )
    return pm_ntv(v).SIMPLICIAL::Bool
end
export is_orbifold


@doc Markdown.doc"""
    issimplicial( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` is simplicial. Hence, this function works just as `is_orbifold`. It is implemented for user convenience.

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
    is_gorenstein( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` is Gorenstein.

# Examples
```jldoctest
julia> is_gorenstein( projective_space( 2 ) )
true
```
"""
function is_gorenstein( v::AbstractNormalToricVariety )
    return pm_ntv(v).GORENSTEIN::Bool
end
export is_gorenstein


@doc Markdown.doc"""
    is_q_gorenstein( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` is Q-Gorenstein.

# Examples
```jldoctest
julia> is_q_gorenstein( projective_space( 2 ) )
true
```
"""
function is_q_gorenstein( v::AbstractNormalToricVariety )
    return pm_ntv(v).Q_GORENSTEIN::Bool
end
export is_q_gorenstein


@doc Markdown.doc"""
    is_fano( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` is fano.

# Examples
```jldoctest
julia> is_fano( projective_space( 2 ) )
true
```
"""
function is_fano( v::AbstractNormalToricVariety )
    return pm_ntv(v).FANO::Bool
end
export is_fano


###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, ntv::AbstractNormalToricVariety)
    # fan = get_polyhedral_fan(ntv)
    pmntv = pm_ntv(ntv)
    ambdim = pmntv.FAN_AMBIENT_DIM
    print(io, "A normal toric variety corresponding to a polyhedral fan in ambient dimension $(ambdim)")
end

