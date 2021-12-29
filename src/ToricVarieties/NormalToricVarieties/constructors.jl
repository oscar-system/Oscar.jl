######################
# 1: The Julia type for ToricVarieties
######################
abstract type AbstractNormalToricVariety end

@attributes mutable struct NormalToricVariety <: AbstractNormalToricVariety
           polymakeNTV::Polymake.BigObject
end
export NormalToricVariety

@attributes mutable struct AffineNormalToricVariety <: AbstractNormalToricVariety
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
A normal, affine, non-complete toric variety
```
"""
function AffineNormalToricVariety(C::Cone)
    # construct the variety
    fan = PolyhedralFan(C)
    pmntv = Polymake.fulton.NormalToricVariety(Oscar.pm_object(fan))
    variety = AffineNormalToricVariety(pmntv, Dict())
    
    # set known attributes
    set_attribute!(variety, :cone, C)
    set_attribute!(variety, :fan, fan)
    set_attribute!(variety, :isaffine, true)
    set_attribute!(variety, :iscomplete, false)
    set_attribute!(variety, :isprojective, false)
    set_attribute!(variety, :isprojective_space, false)
    
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
A normal, affine, non-complete toric variety
```
"""
function NormalToricVariety(C::Cone)
    # construct the variety
    fan = PolyhedralFan(C)
    pmntv = Polymake.fulton.NormalToricVariety(Oscar.pm_object(fan))
    variety = NormalToricVariety(pmntv, Dict())
    
    # set known attributes
    set_attribute!(variety, :fan, fan)
    set_attribute!(variety, :isaffine, true)
    set_attribute!(variety, :iscomplete, false)
    set_attribute!(variety, :isprojective, false)
    set_attribute!(variety, :isprojective_space, false)
    
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
    variety = NormalToricVariety(pmntv, Dict())
    
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
A normal, affine, non-complete toric variety

julia> affineVariety = AffineNormalToricVariety(v)
A normal, affine, non-complete toric variety
```
"""
function AffineNormalToricVariety(v::NormalToricVariety)
    # check input
    isaffine(v) || error("Cannot construct affine toric variety from non-affine input")
    
    # set variety
    variety = AffineNormalToricVariety(pm_object(v), Dict())
    
    # set properties
    set_attribute!(variety, :isaffine, true)
    set_attribute!(variety, :iscomplete, false)
    set_attribute!(variety, :isprojective, false)
    set_attribute!(variety, :isprojective_space, false)    
    
    # construct the affine variety and copy all cached information from v
    return variety
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
A normal, non-affine, smooth, projective, gorenstein, q-gorenstein, fano, 2-dimensional toric variety without torusfactor
```
"""
function toric_projective_space(d::Int)
    # construct the variety
    f = normal_fan(Oscar.simplex(d))
    pm_object = Polymake.fulton.NormalToricVariety(Oscar.pm_object(f))
    variety = NormalToricVariety(pm_object, Dict())
    
    # set properties
    set_attribute!(variety, :isaffine, false)
    set_attribute!(variety, :isprojective, true)
    set_attribute!(variety, :isprojective_space, true)
    set_attribute!(variety, :issmooth, true)
    set_attribute!(variety, :iscomplete, true)
    set_attribute!(variety, :hastorusfactor, false)
    set_attribute!(variety, :isorbifold, true)
    set_attribute!(variety, :issimplicial, true)
    set_attribute!(variety, :isgorenstein, true)
    set_attribute!(variety, :isq_gorenstein, true)
    set_attribute!(variety, :isfano, true)
    
    # set attributes
    set_attribute!(variety, :dim, d)
    set_attribute!(variety, :dim_of_torusfactor, 0)
    set_attribute!(variety, :euler_characteristic, d+1)
    set_attribute!(variety, :character_lattice, abelian_group([0 for i in 1:d]))
    set_attribute!(variety, :torusinvariant_divisor_group, abelian_group([0 for i in 1:d+1]))
    set_attribute!(variety, :map_from_cartier_divisor_group_to_torus_invariant_divisor_group, Hecke.identity_map(torusinvariant_divisor_group(variety)))
    set_attribute!(variety, :map_from_cartier_divisor_group_to_picard_group, map_from_weil_divisors_to_class_group(variety))
    set_attribute!(variety, :stanley_reisner_ideal, ideal([prod(Hecke.gens(cox_ring(variety)))]))
    set_attribute!(variety, :irrelevant_ideal, ideal(Hecke.gens(cox_ring(variety))))
    betti_numbers = [if iseven(i) fmpz(1) else fmpz(0) end for i in 0:2*d]
    set_attribute!(variety, :betti_number, betti_numbers)
    
    # return the variety
    return variety
end
export toric_projective_space


@doc Markdown.doc"""
    hirzebruch_surface(r::Int)

Constructs the r-th Hirzebruch surface.

# Examples
```jldoctest
julia> hirzebruch_surface(5)
A normal, non-affine, smooth, projective, gorenstein, q-gorenstein, non-fano, 2-dimensional toric variety without torusfactor
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
    set_attribute!(variety, :isprojective_space, false)
    set_attribute!(variety, :issmooth, true)
    set_attribute!(variety, :iscomplete, true)
    set_attribute!(variety, :hastorusfactor, false)
    set_attribute!(variety, :isorbifold, true)
    set_attribute!(variety, :issimplicial, true)
    set_attribute!(variety, :isgorenstein, true)
    set_attribute!(variety, :isq_gorenstein, true)
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
    set_attribute!(variety, :map_from_cartier_divisor_group_to_torus_invariant_divisor_group, Hecke.identity_map(torusinvariant_divisor_group(variety)))
    set_attribute!(variety, :map_from_cartier_divisor_group_to_picard_group, map_from_weil_divisors_to_class_group(variety))
    gens = Hecke.gens(cox_ring(variety))
    set_attribute!(variety, :stanley_reisner_ideal, ideal([gens[1]*gens[3],gens[2]*gens[4]]))
    set_attribute!(variety, :irrelevant_ideal, ideal([gens[1]*gens[2], gens[3]*gens[2], gens[1]*gens[4], gens[3]*gens[4]]))
    set_attribute!(variety, :betti_number, [fmpz(1),fmpz(0),fmpz(2),fmpz(0),fmpz(1)])
    
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
A normal, non-affine, smooth, projective, gorenstein, q-gorenstein, fano, 2-dimensional toric variety without torusfactor
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
        return toric_projective_space(2)
    end
    
    # construct the "true" toric del Pezzo surfaces
    if b == 1
        fan_rays = [1 0; 0 1; -1 0; -1 -1]
        cones = IncidenceMatrix([[1,2],[2,3],[3,4],[4,1]])
    end
    if b == 2
        fan_rays = [1 0; 0 1; -1 0; -1 -1; 0 -1]
        cones = IncidenceMatrix([[1,2],[2,3],[3,4],[4,5],[5,1]])
    end
    if b == 3
        fan_rays = [1 0; 1 1; 0 1; -1 0; -1 -1; 0 -1]
        cones = IncidenceMatrix([[1,2],[2,3],[3,4],[4,5],[5,6],[6,1]])
    end
    variety = NormalToricVariety(PolyhedralFan(fan_rays, cones))
    
    # set properties
    set_attribute!(variety, :isaffine, false)
    set_attribute!(variety, :isprojective, true)
    set_attribute!(variety, :isprojective_space, false)
    set_attribute!(variety, :issmooth, true)
    set_attribute!(variety, :iscomplete, true)
    set_attribute!(variety, :hastorusfactor, false)
    set_attribute!(variety, :isorbifold, true)
    set_attribute!(variety, :issimplicial, true)
    set_attribute!(variety, :isgorenstein, true)
    set_attribute!(variety, :isq_gorenstein, true)
    set_attribute!(variety, :isfano, true)
    
    # set attributes that depend on b
    if b == 1
        set_attribute!(variety, :euler_characteristic, 4)
        set_attribute!(variety, :torusinvariant_divisor_group, free_abelian_group(4))    
        ring = PolynomialRing(QQ, :x=>1:2, :e=>1, :x=>3)[1]
        weights = [map_from_weil_divisors_to_class_group(variety)(x) for x in Hecke.gens(torusinvariant_divisor_group(variety))]
        gens = Hecke.gens(cox_ring(variety))
        set_attribute!(variety, :cox_ring, grade(ring,weights)[1])
        set_attribute!(variety, :stanley_reisner_ideal, ideal([gens[1]*gens[3], gens[2]*gens[4]]))
        set_attribute!(variety, :irrelevant_ideal, ideal([gens[3]*gens[4], gens[1]*gens[4], gens[1]*gens[2], gens[3]*gens[2]]))
        set_attribute!(variety, :betti_number, [fmpz(1),fmpz(0),fmpz(2),fmpz(0),fmpz(1)])
    end
    if b == 2
        set_attribute!(variety, :euler_characteristic, 5)
        set_attribute!(variety, :torusinvariant_divisor_group, free_abelian_group(5))
        ring = PolynomialRing(QQ, :x=>1:2, :e=>1, :x=>3, :e=>2)[1]
        weights = [map_from_weil_divisors_to_class_group(variety)(x) for x in Hecke.gens(torusinvariant_divisor_group(variety))]
        gens = Hecke.gens(cox_ring(variety))
        set_attribute!(variety, :cox_ring, grade(ring,weights)[1])
        set_attribute!(variety, :stanley_reisner_ideal, ideal([gens[1]*gens[3], gens[1]*gens[4], 
                                                                                          gens[2]*gens[4], gens[2]*gens[5], gens[3]*gens[5]]))
        set_attribute!(variety, :irrelevant_ideal, ideal([gens[3]*gens[4]*gens[5], gens[1]*gens[4]*gens[5], 
                                                                                 gens[1]*gens[2]*gens[5], gens[1]*gens[2]*gens[3], gens[2]*gens[3]*gens[4]]))
        set_attribute!(variety, :betti_number, [fmpz(1),fmpz(0),fmpz(3),fmpz(0),fmpz(1)])
    end
    if b == 3
        set_attribute!(variety, :euler_characteristic, 6)
        set_attribute!(variety, :torusinvariant_divisor_group, free_abelian_group(6))
        ring = PolynomialRing(QQ, :x=>1, :e=>3, :x=>2, :e=>1, :x=>3, :e=>2)[1]
        weights = [map_from_weil_divisors_to_class_group(variety)(x) for x in Hecke.gens(torusinvariant_divisor_group(variety))]
        gens = Hecke.gens(cox_ring(variety))
        set_attribute!(variety, :cox_ring, grade(ring,weights)[1])
        set_attribute!(variety, :stanley_reisner_ideal, ideal([gens[1]*gens[3], gens[1]*gens[4], gens[1]*gens[5], gens[2]*gens[4], 
                                                                                          gens[2]*gens[5], gens[2]*gens[6], gens[3]*gens[5], gens[3]*gens[6], gens[4]*gens[6]]))
        set_attribute!(variety, :irrelevant_ideal, ideal([gens[3]*gens[4] *gens[5]*gens[6], gens[1]*gens[4] *gens[5]*gens[6],
                                                                                 gens[1]*gens[2] *gens[5]*gens[6], gens[1]*gens[2] *gens[3]*gens[6],
                                                                                 gens[1]*gens[2] *gens[3]*gens[4], gens[2]*gens[3] *gens[4]*gens[5]]))
        set_attribute!(variety, :betti_number, [fmpz(1),fmpz(0),fmpz(4),fmpz(0),fmpz(1)])
    end
    
    # set further attributes
    set_attribute!(variety, :dim, 2)
    set_attribute!(variety, :dim_of_torusfactor, 0)
    set_attribute!(variety, :character_lattice, free_abelian_group(2))
    set_attribute!(variety, :map_from_cartier_divisor_group_to_torus_invariant_divisor_group, Hecke.identity_map(torusinvariant_divisor_group(variety)))
    set_attribute!(variety, :map_from_cartier_divisor_group_to_picard_group, map_from_weil_divisors_to_class_group(variety))
    
    # return the result
    return variety    
end
export del_pezzo


############################
# 4: Advanced constructions
############################

@doc Markdown.doc"""
    blowup_on_ith_minimal_torus_orbit(v::AbstractNormalToricVariety, n::Int)

Return the blowup of the normal toric variety `v` on its i-th minimal torus orbit.

# Examples
```jldoctest
julia> P2 = toric_projective_space(2)
A normal, non-affine, smooth, projective, gorenstein, q-gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> blowup_on_ith_minimal_torus_orbit(P2,1)
A normal toric variety
```
"""
function blowup_on_ith_minimal_torus_orbit(v::AbstractNormalToricVariety, n::Int)
    return NormalToricVariety(starsubdivision(fan(v), n))
end
export blowup_on_ith_minimal_torus_orbit


@doc Markdown.doc"""
    Base.:*(v::AbstractNormalToricVariety, w::AbstractNormalToricVariety)

Return the Cartesian/direct product of two normal toric varieties `v` and `w`.

# Examples
```jldoctest
julia> P2 = toric_projective_space(2)
A normal, non-affine, smooth, projective, gorenstein, q-gorenstein, fano, 2-dimensional toric variety without torusfactor

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

    # affine?
    if has_attribute(v, :isaffine)
        if get_attribute(v, :isaffine)
            push!(properties_string, "affine")
        else
            push!(properties_string, "non-affine")
        end
    end
    
    # smooth/simplicial?
    if has_attribute(v, :issmooth)
        if get_attribute(v, :issmooth)
            push!(properties_string, "smooth")
        else
            if has_attribute(v, :issimplicial)
                if get_attribute(v, :issimplicial)
                    push!(properties_string, "simplicial")
                else
                    push!(properties_string, "non-smooth")
                end
            end
        end
    end
    
    # complete/projective?
    if has_attribute(v, :iscomplete)
        if get_attribute(v, :iscomplete)
            if has_attribute(v, :isprojective)
                if get_attribute(v, :isprojective)
                    push!(properties_string, "projective")
                end
            else
                push!(properties_string, "complete")
            end
        else
            push!(properties_string, "non-complete")
        end
    end
    
    # gorenstein?
    if has_attribute(v, :isgorenstein)
        if get_attribute(v, :isgorenstein)
            push!(properties_string, "gorenstein")
        else
            push!(properties_string, "non-gorenstein")
        end
    end
    
    # q-gorenstein?
    if has_attribute(v, :isq_gorenstein)
        if get_attribute(v, :isq_gorenstein)
            push!(properties_string, "q-gorenstein")
        else
            push!(properties_string, "non-q-gorenstein")
        end
    end
    
    # fano?
    if has_attribute(v, :isfano)
        if get_attribute(v, :isfano)
            push!(properties_string, "fano")
        else
            push!(properties_string, "non-fano")
        end
    end
    
    # dimension?
    if has_attribute(v, :dim)
        push!(properties_string, string(dim(v))*"-dimensional")
    end
    
    # torusfactor?
    if has_attribute(v, :hastorusfactor)
        if get_attribute(v, :hastorusfactor)
            push!(properties_string, "toric variety with torusfactor")
        else
            push!(properties_string, "toric variety without torusfactor")
        end
    else
        push!(properties_string, "toric variety")
    end
    
    # print the information
    join(io, properties_string, ", ", " ")
end
