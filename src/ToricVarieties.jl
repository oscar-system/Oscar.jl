######################
# 1: The Julia type for ToricVarieties
######################
abstract type AbstractNormalToricVariety end

struct NormalToricVariety <: AbstractNormalToricVariety
           GapNTV::GapObj
           polymakeNTV::Polymake.BigObject
end
export NormalToricVariety

struct AffineNormalToricVariety <: AbstractNormalToricVariety
           GapNTV::GapObj
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
    pmntv = Polymake.fulton.NormalToricVariety(Oscar.pm_fan(PF))
    return NormalToricVariety(ntv_polymake2gap(pmntv), pmntv)
end



@doc Markdown.doc"""
    NormalToricVariety( r::Matrix{Int}, c::Vector{Vector{Int}} )

Construct the normal toric variety whose fan has ray generators `r` and maximal cones `c`.

# Examples
```julia-repl
julia> NormalToricVariety( [-1 5; 0 1; 1 0; 0 -1], [[1,2],[2,3],[3,4],[4,1]] )
NormalToricVariety(GAP: <A toric variety of dimension 2>, Polymake.BigObjectAllocated(Ptr{Nothing} @0x0000561365a4f3d0))
```
"""
function NormalToricVariety( rays::Matrix{Int}, cones::Vector{Vector{Int}} )
    # construct the toric variety in GAP
    gap_rays = GapObj( rays, recursive = true )
    gap_cones = GapObj( cones, recursive = true )
    fan = GAP.Globals.Fan( gap_rays, gap_cones )
    variety = GAP.Globals.ToricVariety( fan )

    Incidence = Oscar.IncidenceMatrix(cones)
    arr = @Polymake.convert_to Array{Set{Int}} Polymake.common.rows(Incidence.pm_incidencematrix)

    pmntv = Polymake.fulton.NormalToricVariety(
        RAYS = Oscar.matrix_for_polymake(rays),
        MAXIMAL_CONES = arr,
    )

    # wrap it into a struct and return
    return NormalToricVariety( variety, pmntv )
end

@doc Markdown.doc"""
    NormalToricVariety( v::GapObj )

Construct the Julia wrapper for a `GAP` toric variety `v`.
"""
function NormalToricVariety(GapNTV::GapObj)
   pmNTV = ntv_gap2polymake(GapNTV)
   return NormalToricVariety(GapNTV, pmNTV)
end
export NormalToricVariety

@doc Markdown.doc"""
    ntv_gap2polymake( v::GapObj )

Convert a `GAP` toric variety `v` into a `Polymake` toric variety.
"""
function ntv_gap2polymake(GapNTV::GapObj)
    ff = GAP.Globals.Fan(GapNTV)
    R = GAP.Globals.RayGenerators(ff)
    MC = GAP.Globals.RaysInMaximalCones(ff)
    rays = Matrix{Int}(R)
    cones = [findall(x->x!=0, vec) for vec in [[e for e in mc] for mc in MC]]
    Incidence = Oscar.IncidenceMatrix(cones)
    arr = @Polymake.convert_to Array{Set{Int}} Polymake.common.rows(Incidence.pm_incidencematrix)
    pmntv = Polymake.fulton.NormalToricVariety(
        RAYS = Oscar.matrix_for_polymake(rays),
        MAXIMAL_CONES = arr,
    )
    return pmntv
end
export ntv_gap2polymake

@doc Markdown.doc"""
    ntv_polymake2gap( v::Polymake.BigObject )

Convert a `Polymake` toric variety `v` into a `GAP` toric variety.
"""
function ntv_polymake2gap(polymakeNTV::Polymake.BigObject)
    rays = Matrix{Int}(polymakeNTV.RAYS)
    gap_rays = GapObj( rays, recursive = true )
    cones = [findall(x->x!=0, v) for v in eachrow(polymakeNTV.MAXIMAL_CONES)]
    gap_cones = GapObj( cones, recursive = true )
    fan = Oscar.GAP.Globals.Fan( gap_rays, gap_cones )
    variety = Oscar.GAP.Globals.ToricVariety( fan )
    return variety
end
export ntv_polymake2gap


######################
# 3: Special constructors
######################

@doc Markdown.doc"""
    projective_space( d::Int )

Construct the projective space of dimension `d`.

# Examples
```julia-repl
julia> projective_space( 2 )
NormalToricVariety(GAP: <A projective toric variety of dimension 2>, Polymake.BigObjectAllocated(Ptr{Nothing} @0x0000562fbcc21b70))
```
"""
function projective_space( d::Int )
    # construct the projective space in gap
    variety = GAP.Globals.ProjectiveSpace( d )
    f = Polymake.fan.normal_fan(Polymake.polytope.simplex(d))
    pmntv = Polymake.fulton.NormalToricVariety(f)
    
    # wrap it and return
    return NormalToricVariety( variety, pmntv )
end
export projective_space


@doc Markdown.doc"""
    hirzebruch_surface( r::Int )

Constructs the r-th Hirzebruch surface.

# Examples
```julia-repl
julia> hirzebruch_surface( 5 )
NormalToricVariety(GAP: <A toric variety of dimension 2>, Polymake.BigObjectAllocated(Ptr{Nothing} @0x0000562fbcc21b70))
```
"""
function hirzebruch_surface( r::Int )
    pmntv = Polymake.fulton.hirzebruch_surface(r)
    return NormalToricVariety(ntv_polymake2gap(pmntv), pmntv)
end
export hirzebruch_surface


@doc Markdown.doc"""
    delPezzo( b::Int )

Constructs the delPezzo surface with b blowups for b at most 3.

# Examples
```julia-repl
julia> del_pezzo( 3 )
NormalToricVariety(GAP: <A toric variety of dimension 2>, Polymake.BigObjectAllocated(Ptr{Nothing} @0x0000562fbcc21b70))
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
```julia-repl
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
```julia-repl
julia> isaffine( projective_space( 2 ) )
false
```
"""
function isaffine( v::AbstractNormalToricVariety )
    return pm_ntv(v).AFFINE
end
export isaffine


@doc Markdown.doc"""
    isprojective( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` is projective, i.e. if the fan of `v` is the the normal fan of a polytope.

# Examples
```julia-repl
julia> isprojective( projective_space( 2 ) )
true
```
"""
function isprojective( v::AbstractNormalToricVariety )
    return pm_ntv(v).PROJECTIVE
end
export isprojective


@doc Markdown.doc"""
    issmooth( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` is smooth.

# Examples
```julia-repl
julia> issmooth( projective_space( 2 ) )
true
```
"""
function issmooth( v::AbstractNormalToricVariety )
    return pm_ntv(v).SMOOTH
end
export issmooth


@doc Markdown.doc"""
    iscomplete( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` is complete.

# Examples
```julia-repl
julia> iscomplete( projective_space( 2 ) )
true
```
"""
function iscomplete( v::AbstractNormalToricVariety )
    return pm_ntv(v).COMPLETE
end
export iscomplete


@doc Markdown.doc"""
    has_torusfactor( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` has a torus factor.

# Examples
```julia-repl
julia> has_torusfactor( projective_space( 2 ) )
false
```
"""
function has_torusfactor( v::AbstractNormalToricVariety )
    return GAP.Globals.HasTorusfactor( v.GapNTV )::Bool
end
export has_torusfactor


@doc Markdown.doc"""
    is_orbifold( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` is an orbifold.

# Examples
```julia-repl
julia> is_orbifold( projective_space( 2 ) )
true
```
"""
function is_orbifold( v::AbstractNormalToricVariety )
    return GAP.Globals.IsOrbifold( v.GapNTV )::Bool
end
export is_orbifold


@doc Markdown.doc"""
    issimplicial( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` is simplicial. Hence, this function works just as `is_orbifold`. It is implemented for user convenience.

# Examples
```julia-repl
julia> issimplicial( projective_space( 2 ) )
true
```
"""
function issimplicial( v::AbstractNormalToricVariety )
    return pm_ntv(v).SIMPLICIAL
end
export issimplicial


@doc Markdown.doc"""
    is_isomorphic_to_projective_space( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` is isomorphic to projective space.

# Examples
```julia-repl
julia> is_isomorphic_to_projective_space( projective_space( 2 ) )
true
```
"""
function is_isomorphic_to_projective_space( v::AbstractNormalToricVariety )
    return GAP.Globals.IsIsomorphicToProjectiveSpace( v.GapNTV )::Bool
end
export is_isomorphic_to_projective_space


@doc Markdown.doc"""
    is_direct_product_of_projective_spaces( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` is isomorphic to a direct product of projective space.

# Examples
```julia-repl
julia> is_direct_product_of_projective_spaces( projective_space( 2 ) )
true
```
"""
function is_direct_product_of_projective_spaces( v::AbstractNormalToricVariety )
    return GAP.Globals.IsDirectProductOfPNs( v.GapNTV )::Bool
end
export is_direct_product_of_projective_spaces


@doc Markdown.doc"""
    is_gorenstein( v::AbstractNormalToricVariety )

Checks if the normal toric variety `v` is Gorenstein.

# Examples
```julia-repl
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
```julia-repl
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
```julia-repl
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

