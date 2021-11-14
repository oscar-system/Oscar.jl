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

function pm_object(v::AbstractNormalToricVariety)
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
julia> C = positive_hull([1 0; 0 1])
A polyhedral cone in ambient dimension 2

julia> antv = AffineNormalToricVariety(C)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""
function AffineNormalToricVariety(C::Cone)
    fan = PolyhedralFan(C)
    pmntv = Polymake.fulton.NormalToricVariety(Oscar.pm_object(fan))
    return AffineNormalToricVariety(pmntv)
end


@doc Markdown.doc"""
    NormalToricVariety(C::Cone)

Construct the (affine) normal toric variety $X_{\Sigma}$ corresponding to a
polyhedral fan $\Sigma = C$ consisting only of the cone `C`.

# Examples
Set `C` to be the positive orthant in two dimensions.
```jldoctest
julia> C = positive_hull([1 0; 0 1])
A polyhedral cone in ambient dimension 2
julia> ntv = NormalToricVariety(C)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""
function NormalToricVariety(C::Cone)
    fan = PolyhedralFan(C)
    pmntv = Polymake.fulton.NormalToricVariety(Oscar.pm_object(fan))
    return NormalToricVariety(pmntv)
end


@doc Markdown.doc"""
    NormalToricVariety(PF::PolyhedralFan)

Construct the normal toric variety $X_{PF}$ corresponding to a polyhedral fan `PF`.

# Examples
Take `PF` to be the normal fan of the square.
```jldoctest
julia> square = cube(2)
A polyhedron in ambient dimension 2

julia> nf = normal_fan(square)
A polyhedral fan in ambient dimension 2

julia> ntv = NormalToricVariety(nf)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""    
function NormalToricVariety(PF::PolyhedralFan)
    fan = Oscar.pm_object(PF)
    pmntv = Polymake.fulton.NormalToricVariety(fan)
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
julia> square = cube(2)
A polyhedron in ambient dimension 2

julia> ntv = NormalToricVariety(square)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""    
function NormalToricVariety(P::Polyhedron)
    fan = normal_fan(P)
    return NormalToricVariety(fan)
end

export NormalToricVariety

@doc Markdown.doc"""
    AffineNormalToricVariety(v::NormalToricVariety)

For internal design, we make a strict distinction between
normal toric varieties and affine toric varieties.
Given an affine, normal toric variety `v`,
this method turns it into an affine toric variety.

# Examples
```jldoctest
julia> v = NormalToricVariety(positive_hull([1 0; 0 1]))
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> affineVariety = AffineNormalToricVariety(v)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""
function AffineNormalToricVariety(v::NormalToricVariety)
    isaffine(v) || error("Cannot construct affine toric variety from non-affine input")
    return AffineNormalToricVariety(pm_object(v))
end



######################
# 3: Special constructors
######################

@doc Markdown.doc"""
    toric_projective_space(d::Int)

Construct the projective space of dimension `d`.

# Examples
```jldoctest
julia> toric_projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""
function toric_projective_space(d::Int)
    f = normal_fan(Oscar.simplex(d))
    pm_object = Polymake.fulton.NormalToricVariety(Oscar.pm_object(f))
    return NormalToricVariety(pm_object)
end
export toric_projective_space


@doc Markdown.doc"""
    hirzebruch_surface(r::Int)

Constructs the r-th Hirzebruch surface.

# Examples
```jldoctest
julia> hirzebruch_surface(5)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""
function hirzebruch_surface(r::Int)
    rays = [1 0; 0 1; -1 r; 0 -1]
    cones = IncidenceMatrix([[1,2],[2,3],[3,4],[4,1]])
    return NormalToricVariety(PolyhedralFan(rays, cones))
end
export hirzebruch_surface


@doc Markdown.doc"""
    del_pezzo(b::Int)

Constructs the delPezzo surface with b blowups for b at most 3.

# Examples
```jldoctest
julia> del_pezzo(3)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""
function del_pezzo(b::Int)
    if b < 0
        throw(ArgumentError("Number of blowups for construction of delPezzo surfaces must be non-negative."))
    end
    if b == 0 
        return toric_projective_space(2)
    end
    if b == 1
        rays = [1 0; 0 1; -1 0; -1 -1]
        cones = IncidenceMatrix([[1,2],[2,3],[3,4],[4,1]])
        return NormalToricVariety(PolyhedralFan(rays, cones))
    end
    if b == 2
        rays = [1 0; 0 1; -1 0; -1 -1; 0 -1]
        cones = IncidenceMatrix([[1,2],[2,3],[3,4],[4,5],[5,1]])
        return NormalToricVariety(PolyhedralFan(rays, cones))
    end
    if b == 3
        rays = [1 0; 1 1; 0 1; -1 0; -1 -1; 0 -1]
        cones = IncidenceMatrix([[1,2],[2,3],[3,4],[4,5],[5,6],[6,1]])
        return NormalToricVariety(PolyhedralFan(rays, cones))
    end
    if b > 3
        throw(ArgumentError("delPezzo surfaces with more than 3 blowups are realized as subvarieties of toric ambient spaces. This is currently not supported."))
    end
end
export del_pezzo


############################
# 4: Advanced constructions
############################

@doc Markdown.doc"""
    blowup_on_ith_minimal_torus_orbit(v::AbstractNormalToricVariety, n::Int)

Computes the blowup of the normal toric variety `v` on its i-th minimal torus orbit.

# Examples
```jldoctest
julia> P2 = toric_projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> blowup_on_ith_minimal_torus_orbit(P2,1)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""
function blowup_on_ith_minimal_torus_orbit(v::AbstractNormalToricVariety, n::Int)
    return NormalToricVariety( starsubdivision( fan_of_variety( v ), n ) )
end
export blowup_on_ith_minimal_torus_orbit


@doc Markdown.doc"""
    Base.:*(v::AbstractNormalToricVariety, w::AbstractNormalToricVariety)

Computes the Cartesian/direct product of two normal toric varieties `v` and `w`.

# Examples
```jldoctest
julia> P2 = toric_projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> P2 * P2
A normal toric variety corresponding to a polyhedral fan in ambient dimension 4
```
"""
function Base.:*(v::AbstractNormalToricVariety, w::AbstractNormalToricVariety)
    return NormalToricVariety(fan_of_variety(v)*fan_of_variety(w))
end


###############################################################################
###############################################################################
### 5: Display
###############################################################################
###############################################################################
function Base.show(io::IO, ntv::AbstractNormalToricVariety)
    # fan = get_polyhedral_fan(ntv)
    pmntv = pm_object(ntv)
    ambdim = pmntv.FAN_AMBIENT_DIM
    print(io, "A normal toric variety corresponding to a polyhedral fan in ambient dimension $(ambdim)")
end
