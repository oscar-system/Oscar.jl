import Oscar.projective_space

######################
# 1: The Julia type for ToricVarieties
######################
abstract type AbstractNormalToricVariety <: _FanLikeType end

@attributes mutable struct NormalToricVariety <: AbstractNormalToricVariety
           polymakeNTV::Polymake.BigObject
           NormalToricVariety(polymakeNTV::Polymake.BigObject) = new(polymakeNTV)
end
export NormalToricVariety

@attributes mutable struct AffineNormalToricVariety <: AbstractNormalToricVariety
           polymakeNTV::Polymake.BigObject
           AffineNormalToricVariety(polymakeNTV::Polymake.BigObject) = new(polymakeNTV)
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
A normal, affine toric variety
```
"""
function AffineNormalToricVariety(C::Cone)
    # construct the variety
    fan = PolyhedralFan(C)
    pmntv = Polymake.fulton.NormalToricVariety(Oscar.pm_object(fan))
    variety = AffineNormalToricVariety(pmntv)
    
    # set known attributes
    set_attribute!(variety, :cone, C)
    set_attribute!(variety, :fan, fan)
    set_attribute!(variety, :isaffine, true)
    set_attribute!(variety, :iscomplete, false)
    set_attribute!(variety, :isprojective, false)
    set_attribute!(variety, :is_projective_space, false)
    set_attribute!(variety, :picard_group, free_abelian_group(0))
    
    # return
    return variety
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
A normal, affine toric variety
```
"""
function NormalToricVariety(C::Cone)
    # construct the variety
    fan = PolyhedralFan(C)
    pmntv = Polymake.fulton.NormalToricVariety(Oscar.pm_object(fan))
    variety = NormalToricVariety(pmntv)
    
    # set known attributes
    set_attribute!(variety, :fan, fan)
    set_attribute!(variety, :isaffine, true)
    set_attribute!(variety, :iscomplete, false)
    set_attribute!(variety, :isprojective, false)
    set_attribute!(variety, :is_projective_space, false)
    set_attribute!(variety, :picard_group, free_abelian_group(0))
    
    # return
    return variety
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
A normal toric variety
```
"""    
function NormalToricVariety(PF::PolyhedralFan)
    # construct the variety    
    fan = Oscar.pm_object(PF)
    pmntv = Polymake.fulton.NormalToricVariety(fan)
    variety = NormalToricVariety(pmntv)
    
    # set attributes
    set_attribute!(variety, :fan, PF)
    
    # return
    return variety
end


@doc Markdown.doc"""
    NormalToricVariety(P::Polyhedron)

Construct the normal toric variety $X_{\Sigma_P}$ corresponding to the normal
fan $\Sigma_P$ of the given polyhedron `P`.

Note that this only coincides with the projective variety associated to `P`
from the affine relations of the lattice points in `P`, if `P` is very ample.

# Examples
Set `P` to be a square.
```jldoctest
julia> square = cube(2)
A polyhedron in ambient dimension 2

julia> ntv = NormalToricVariety(square)
A normal toric variety
```
"""    
function NormalToricVariety(P::Polyhedron)
    fan = normal_fan(P)
    variety = NormalToricVariety(fan)
    set_attribute!(variety, :polyhedron, P)
    return variety
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
A normal, affine toric variety

julia> affineVariety = AffineNormalToricVariety(v)
A normal, affine toric variety
```
"""
function AffineNormalToricVariety(v::NormalToricVariety)
    # check input
    isaffine(v) || error("Cannot construct affine toric variety from non-affine input")
    
    # set variety
    variety = AffineNormalToricVariety(pm_object(v))
    
    # set properties of variety
    set_attribute!(variety, :isaffine, true)
    set_attribute!(variety, :iscomplete, false)
    set_attribute!(variety, :isprojective, false)
    set_attribute!(variety, :is_projective_space, false)
    
    # construct the affine variety and copy all cached information from v
    return variety
end


######################
# 3: Special constructors
######################

@doc Markdown.doc"""
    affine_space(::Type{NormalToricVariety}, d::Int)

Constructs the (toric) affine space of dimension `d`.

# Examples
```jldoctest
julia> affine_space(NormalToricVariety, 2)
A normal, affine, 2-dimensional toric variety
```
"""
function affine_space(::Type{NormalToricVariety}, d::Int)
    # construct the cone of the variety
    m = zeros(Int, d, d)
    for i in 1:d
        m[i,i] = 1
    end
    C = positive_hull(m)
    
    # construct the variety
    fan = PolyhedralFan(C)
    pmntv = Polymake.fulton.NormalToricVariety(Oscar.pm_object(fan))
    variety = NormalToricVariety(pmntv)
    
    # set known properties
    set_attribute!(variety, :isaffine, true)
    set_attribute!(variety, :iscomplete, false)
    set_attribute!(variety, :isprojective, false)
    set_attribute!(variety, :isprojective_space, false)
    
    # set attributes
    set_attribute!(variety, :fan, fan)
    set_attribute!(variety, :dim, d)
    set_attribute!(variety, :dim_of_torusfactor, 0)
    
    # return the variety
    return variety
end
export affine_space


@doc Markdown.doc"""
    projective_space(::Type{NormalToricVariety}, d::Int)

Construct the projective space of dimension `d`.

# Examples
```jldoctest
julia> projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor
```
"""
function projective_space(::Type{NormalToricVariety}, d::Int)
    # construct the variety
    f = normal_fan(Oscar.simplex(d))
    pm_object = Polymake.fulton.NormalToricVariety(Oscar.pm_object(f))
    variety = NormalToricVariety(pm_object)
    
    # set properties
    set_attribute!(variety, :isaffine, false)
    set_attribute!(variety, :isprojective, true)
    set_attribute!(variety, :is_projective_space, true)
    set_attribute!(variety, :issmooth, true)
    set_attribute!(variety, :iscomplete, true)
    set_attribute!(variety, :hastorusfactor, false)
    set_attribute!(variety, :isorbifold, true)
    set_attribute!(variety, :issimplicial, true)
    set_attribute!(variety, :isgorenstein, true)
    set_attribute!(variety, :is_q_gorenstein, true)
    set_attribute!(variety, :isfano, true)
    
    # set attributes
    set_attribute!(variety, :dim, d)
    set_attribute!(variety, :dim_of_torusfactor, 0)
    set_attribute!(variety, :euler_characteristic, d+1)
    set_attribute!(variety, :character_lattice, free_abelian_group(d))
    set_attribute!(variety, :torusinvariant_divisor_group, free_abelian_group(d+1))
    set_attribute!(variety, :map_from_torus_invariant_cartier_divisor_group_to_torus_invariant_weil_divisor_group, identity_map(torusinvariant_divisor_group(variety)))
    set_attribute!(variety, :map_from_torus_invariant_cartier_divisor_group_to_picard_group, map_from_torus_invariant_weil_divisor_group_to_class_group(variety))
    betti_numbers = fill(fmpz(1), d+1)
    set_attribute!(variety, :betti_number, betti_numbers)
    
    # return the variety
    return variety
end
export projective_space


@doc Markdown.doc"""
    hirzebruch_surface(r::Int)

Constructs the r-th Hirzebruch surface.

# Examples
```jldoctest
julia> hirzebruch_surface(5)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor
```
"""
function hirzebruch_surface(r::Int)
    # construct the variety
    fan_rays = [1 0; 0 1; -1 r; 0 -1]
    cones = IncidenceMatrix([[1,2],[2,3],[3,4],[4,1]])
    variety = NormalToricVariety(PolyhedralFan(fan_rays, cones))
    
    # set properties
    set_attribute!(variety, :isaffine, false)
    set_attribute!(variety, :isprojective, true)
    set_attribute!(variety, :is_projective_space, false)
    set_attribute!(variety, :issmooth, true)
    set_attribute!(variety, :iscomplete, true)
    set_attribute!(variety, :hastorusfactor, false)
    set_attribute!(variety, :isorbifold, true)
    set_attribute!(variety, :issimplicial, true)
    set_attribute!(variety, :isgorenstein, true)
    set_attribute!(variety, :is_q_gorenstein, true)
    if abs(r) <= 1
        set_attribute!(variety, :isfano, true)
    else
        set_attribute!(variety, :isfano, false)
    end
    
    # set attributes
    set_attribute!(variety, :dim, 2)
    set_attribute!(variety, :dim_of_torusfactor, 0)
    set_attribute!(variety, :euler_characteristic, 4)
    set_attribute!(variety, :character_lattice, free_abelian_group(2))
    set_attribute!(variety, :torusinvariant_divisor_group, free_abelian_group(4))
    set_attribute!(variety, :map_from_torus_invariant_cartier_divisor_group_to_torus_invariant_weil_divisor_group, identity_map(torusinvariant_divisor_group(variety)))
    set_attribute!(variety, :map_from_torus_invariant_cartier_divisor_group_to_picard_group, map_from_torus_invariant_weil_divisor_group_to_class_group(variety))
    set_attribute!(variety, :betti_number, [fmpz(1),fmpz(2),fmpz(1)])
    
    # return the result
    return variety
end
export hirzebruch_surface


@doc Markdown.doc"""
    del_pezzo(b::Int)

Constructs the delPezzo surface with b blowups for b at most 3.

# Examples
```jldoctest
julia> del_pezzo(3)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor
```
"""
function del_pezzo(b::Int)
    # check for valid input
    if b < 0
        throw(ArgumentError("Number of blowups for construction of delPezzo surfaces must be non-negative."))
    end
    if b > 3
        throw(ArgumentError("delPezzo surfaces with more than 3 blowups are realized as subvarieties of toric ambient spaces. This is currently not supported."))
    end
    
    # special case of projective space
    if b == 0 
        return projective_space(NormalToricVariety, 2)
    end
    
    # construct the "true" toric del Pezzo surfaces
    if b == 1
        fan_rays = [1 0; 0 1; -1 -1; -1 0]
        cones = IncidenceMatrix([[1,2],[2,4],[4,3],[3,1]])
    end
    if b == 2
        fan_rays = [1 0; 0 1; -1 -1; -1 0; 0 -1]
        cones = IncidenceMatrix([[1,2],[2,4],[4,3],[3,5],[5,1]])
    end
    if b == 3
        fan_rays = [1 0; 0 1; -1 -1; -1 0; 0 -1; 1 1]
        cones = IncidenceMatrix([[1,6],[6,2],[2,4],[4,3],[3,5],[5,1]])
    end
    variety = NormalToricVariety(PolyhedralFan(fan_rays, cones))
    
    # set properties
    set_attribute!(variety, :isaffine, false)
    set_attribute!(variety, :isprojective, true)
    set_attribute!(variety, :is_projective_space, false)
    set_attribute!(variety, :issmooth, true)
    set_attribute!(variety, :iscomplete, true)
    set_attribute!(variety, :hastorusfactor, false)
    set_attribute!(variety, :isorbifold, true)
    set_attribute!(variety, :issimplicial, true)
    set_attribute!(variety, :isgorenstein, true)
    set_attribute!(variety, :is_q_gorenstein, true)
    set_attribute!(variety, :isfano, true)
    
    # set attributes that depend on b
    if b == 1
        set_attribute!(variety, :euler_characteristic, 4)
        set_attribute!(variety, :torusinvariant_divisor_group, free_abelian_group(4))
        set_attribute!(variety, :coordinate_names, ["x1", "x2", "x3", "e1"])
        set_attribute!(variety, :betti_number, [fmpz(1),fmpz(2),fmpz(1)])
    end
    if b == 2
        set_attribute!(variety, :euler_characteristic, 5)
        set_attribute!(variety, :torusinvariant_divisor_group, free_abelian_group(5))
        set_attribute!(variety, :coordinate_names, ["x1", "x2", "x3", "e1", "e2"])
        set_attribute!(variety, :betti_number, [fmpz(1),fmpz(3),fmpz(1)])
    end
    if b == 3
        set_attribute!(variety, :euler_characteristic, 6)
        set_attribute!(variety, :torusinvariant_divisor_group, free_abelian_group(6))
        set_attribute!(variety, :coordinate_names, ["x1", "x2", "x3", "e1", "e2", "e3"])
        set_attribute!(variety, :betti_number, [fmpz(1),fmpz(4),fmpz(1)])
    end
    
    # set further attributes
    set_attribute!(variety, :dim, 2)
    set_attribute!(variety, :dim_of_torusfactor, 0)
    set_attribute!(variety, :character_lattice, free_abelian_group(2))
    set_attribute!(variety, :map_from_torus_invariant_cartier_divisor_group_to_torus_invariant_weil_divisor_group, identity_map(torusinvariant_divisor_group(variety)))
    set_attribute!(variety, :map_from_torus_invariant_cartier_divisor_group_to_picard_group, map_from_torus_invariant_weil_divisor_group_to_class_group(variety))
    
    # return the result
    return variety    
end
export del_pezzo


############################
# 4: Advanced constructions
############################

@doc Markdown.doc"""
    blowup_on_ith_minimal_torus_orbit(v::AbstractNormalToricVariety, n::Int, coordinate_name::String)

Return the blowup of the normal toric variety `v` on its i-th minimal torus orbit.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> bP2 = blowup_on_ith_minimal_torus_orbit(P2,1,"e")
A normal toric variety over QQ

julia> cox_ring(bP2)
Multivariate Polynomial Ring in x2, x3, x1, e over Rational Field graded by 
  x2 -> [1 0]
  x3 -> [0 1]
  x1 -> [1 0]
  e -> [-1 1]
```
"""
function blowup_on_ith_minimal_torus_orbit(v::AbstractNormalToricVariety, n::Int, coordinate_name::String)
    # compute the blow-up variety
    new_fan = starsubdivision(fan(v), n)
    new_variety = NormalToricVariety(new_fan)
    
    # extract the old and new rays
    # the new cones are in general given by first (in general) permuting the old rays and then adding a new ray (not necessarily at the last position)
    old_rays = rays(fan(v))
    new_rays = rays(new_fan)
    
    # check for name clash with variable name chosen for blowup
    old_vars = [string(x) for x in gens(cox_ring(v))]
    isnothing(findfirst(x->occursin(coordinate_name, x), old_vars)) ||
        throw(ArgumentError("The provided name for the blowup coordinate is already taken as homogeneous coordinate of the provided toric variety."))
    
    # set up Cox ring of new variety
    new_vars = [if new_rays[i] in old_rays old_vars[findfirst(x->x==new_rays[i], old_rays)] else coordinate_name end for i in 1:length(new_rays)]
    set_attribute!(new_variety, :coordinate_names, new_vars)
    weights = [map_from_torus_invariant_weil_divisor_group_to_class_group(new_variety)(x) for x in gens(torusinvariant_divisor_group(new_variety))]
    set_attribute!(new_variety, :cox_ring_weights, weights)
    if has_attribute(v, :coefficient_ring)
        set_attribute!(new_variety, :coefficient_ring, coefficient_ring(v))
    end
    
    # return variety
    return new_variety
end
export blowup_on_ith_minimal_torus_orbit


@doc Markdown.doc"""
    Base.:*(v::AbstractNormalToricVariety, w::AbstractNormalToricVariety)

Return the Cartesian/direct product of two normal toric varieties `v` and `w`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> P2 * P2
A normal toric variety
```
"""
function Base.:*(v::AbstractNormalToricVariety, w::AbstractNormalToricVariety)
    return NormalToricVariety(fan(v)*fan(w))
end


############################
### 5: Display
############################
function Base.show(io::IO, v::AbstractNormalToricVariety)
    # initiate properties string
    properties_string = ["A normal"]

    affine = push_attribute_if_exists!(properties_string, v, :isaffine, "affine")

    simplicial_cb!(a,b) = push_attribute_if_exists!(a, b, :isorbifold, "simplicial")
    push_attribute_if_exists!(properties_string, v, :issmooth, "smooth"; callback=simplicial_cb!)

    projective = nothing
    if isnothing(affine) || !affine
        complete_cb!(a,b) = push_attribute_if_exists!(a, b, :iscomplete, "complete")
        projective = push_attribute_if_exists!(properties_string, v, :isprojective, "projective"; callback=complete_cb!)
    end

    q_gor_cb!(a,b) = push_attribute_if_exists!(a, b, :is_q_gorenstein, "q-gorenstein")
    gorenstein = push_attribute_if_exists!(properties_string, v, :isgorenstein, "gorenstein"; callback=q_gor_cb!)
    
    push_attribute_if_exists!(properties_string, v, :isfano, "fano")
    
    # dimension?
    if has_attribute(v, :dim)
        push!(properties_string, string(dim(v))*"-dimensional")
    end
    
    # join with ","
    properties_string = [join(properties_string, ", ")]
    
    # add coefficient ring and torusfactor
    if has_attribute(v, :coefficient_ring)
        if coefficient_ring(v) == QQ
            push!(properties_string, "toric variety over QQ")
        elseif coefficient_ring(v) == ZZ
            push!(properties_string, "toric variety over ZZ")
        else
            push!(properties_string, "toric variety over coefficient_ring(variety)")
        end
    else
        push!(properties_string, "toric variety")
    end
    push_attribute_if_exists!(properties_string, v, :hastorusfactor, "with torusfactor", "without torusfactor")
    
    # print string
    join(io, properties_string, " ")
end
