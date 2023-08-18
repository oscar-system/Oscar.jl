######################
# Julia type for ToricVarieties
######################
abstract type AbstractNormalToricVariety end

@attributes mutable struct NormalToricVariety <: AbstractNormalToricVariety
           polymakeNTV::Polymake.BigObject
           NormalToricVariety(polymakeNTV::Polymake.BigObject) = new(polymakeNTV)
end

@attributes mutable struct AffineNormalToricVariety <: AbstractNormalToricVariety
           polymakeNTV::Polymake.BigObject
           AffineNormalToricVariety(polymakeNTV::Polymake.BigObject) = new(polymakeNTV)
end

pm_object(v::AbstractNormalToricVariety) = v.polymakeNTV

coefficient_field(::AbstractNormalToricVariety) = QQ

######################
# Constructors
######################

@doc raw"""
    affine_normal_toric_variety(C::Cone; set_attributes::Bool = true)

Construct the affine normal toric variety $U_{C}$ corresponding to a polyhedral
cone `C`.

# Examples
Set `C` to be the positive orthant in two dimensions.
```jldoctest
julia> C = positive_hull([1 0; 0 1])
Polyhedral cone in ambient dimension 2

julia> antv = affine_normal_toric_variety(C)
Normal, affine toric variety
```
"""
function affine_normal_toric_variety(C::Cone; set_attributes::Bool = true)
    fan = polyhedral_fan(C)
    pmntv = Polymake.fulton.NormalToricVariety(Oscar.pm_object(fan))
    variety = AffineNormalToricVariety(pmntv)
    
    if set_attributes
        set_attribute!(variety, :cone, C)
        set_attribute!(variety, :fan, fan)
        set_attribute!(variety, :is_affine, true)
        set_attribute!(variety, :is_complete, false)
        set_attribute!(variety, :is_projective, false)
        set_attribute!(variety, :is_projective_space, false)
        set_attribute!(variety, :picard_group, free_abelian_group(0))
    end
    
    return variety
end


@doc raw"""
    normal_toric_variety(C::Cone; set_attributes::Bool = true)

Construct the (affine) normal toric variety $X_{\Sigma}$ corresponding to a
polyhedral fan $\Sigma = C$ consisting only of the cone `C`.

# Examples
Set `C` to be the positive orthant in two dimensions.
```jldoctest
julia> C = positive_hull([1 0; 0 1])
Polyhedral cone in ambient dimension 2

julia> ntv = normal_toric_variety(C)
Normal, affine toric variety
```
"""
function normal_toric_variety(C::Cone; set_attributes::Bool = true)
    fan = polyhedral_fan(C)
    pmntv = Polymake.fulton.NormalToricVariety(Oscar.pm_object(fan))
    variety = NormalToricVariety(pmntv)
    
    if set_attributes
        set_attribute!(variety, :is_affine, true)
        set_attribute!(variety, :is_projective_space, false)
        set_attribute!(variety, :picard_group, free_abelian_group(0))
    end
    
    return variety
end


@doc raw"""
    normal_toric_variety(rays::AbstractMatrix, max_cones::Vector{Vector{Int64}})

Construct a normal toric variety $X$ by providing the rays and maximal cones
as vector of vectors. By default, this method assumes that the input is not
non-redundant (e.g. that a ray was entered twice by accident). If the user
is certain that no redundancy exists in the entered information, one can
pass `non_redundant = true` as third argument. This will bypass these consistency
checks. In addition, this will ensure that the order of the rays is not
altered by the constructor.

# Examples
```jldoctest
julia> ray_generators = [[1,0], [0, 1], [-1, 5], [0, -1]]
4-element Vector{Vector{Int64}}:
 [1, 0]
 [0, 1]
 [-1, 5]
 [0, -1]

julia> max_cones = [[1, 2], [2, 3], [3, 4], [4, 1]]
4-element Vector{Vector{Int64}}:
 [1, 2]
 [2, 3]
 [3, 4]
 [4, 1]

julia> normal_toric_variety(ray_generators, max_cones)
Normal toric variety

julia> normal_toric_variety(ray_generators, max_cones; non_redundant = true)
Normal toric variety
```
"""
normal_toric_variety(rays::AbstractCollection[RayVector], max_cones::Vector{Vector{Int64}}; non_redundant::Bool = false) = normal_toric_variety(rays, IncidenceMatrix(max_cones); non_redundant = non_redundant)
function normal_toric_variety(rays::AbstractCollection[RayVector], max_cones::IncidenceMatrix; non_redundant::Bool = false)
  fan = polyhedral_fan(rays, max_cones; non_redundant=non_redundant)
  return normal_toric_variety(fan)
end


@doc raw"""
    normal_toric_variety(PF::PolyhedralFan)

Construct the normal toric variety $X_{PF}$ corresponding to a polyhedral fan `PF`.

# Examples
Take `PF` to be the normal fan of the square.
```jldoctest
julia> square = cube(2)
Polyhedron in ambient dimension 2

julia> nf = normal_fan(square)
Polyhedral fan in ambient dimension 2

julia> ntv = normal_toric_variety(nf)
Normal toric variety
```
"""
function normal_toric_variety(PF::PolyhedralFan)
    fan = Oscar.pm_object(PF)
    pmntv = Polymake.fulton.NormalToricVariety(fan)
    variety = NormalToricVariety(pmntv)
    return variety
end


@doc raw"""
    normal_toric_variety(P::Polyhedron; set_attributes::Bool = true)

Construct the normal toric variety $X_{\Sigma_P}$ corresponding to the normal
fan $\Sigma_P$ of the given polyhedron `P`.

Note that this only coincides with the projective variety associated to `P`
from the affine relations of the lattice points in `P`, if `P` is very ample.

# Examples
Set `P` to be a square.
```jldoctest
julia> square = cube(2)
Polyhedron in ambient dimension 2

julia> ntv = normal_toric_variety(square)
Normal toric variety
```
"""
function normal_toric_variety(P::Polyhedron; set_attributes::Bool = true)
    variety = normal_toric_variety(normal_fan(P))
    if set_attributes
        set_attribute!(variety, :polyhedron, P)
    end
    return variety
end


@doc raw"""
    affine_normal_toric_variety(v::NormalToricVariety; set_attributes::Bool = true)

For internal design, we make a strict distinction between
normal toric varieties and affine toric varieties.
Given an affine, normal toric variety `v`,
this method turns it into an affine toric variety.

# Examples
```jldoctest
julia> v = normal_toric_variety(positive_hull([1 0; 0 1]))
Normal, affine toric variety

julia> affineVariety = affine_normal_toric_variety(v)
Normal, affine toric variety
```
"""
function affine_normal_toric_variety(v::NormalToricVariety; set_attributes::Bool = true)
    is_affine(v) || error("Cannot construct affine toric variety from non-affine input")
    variety = AffineNormalToricVariety(pm_object(v))
    
    if set_attributes
        set_attribute!(variety, :is_affine, true)
        set_attribute!(variety, :is_projective_space, false)
    end
    
    return variety
end



######################
# Display
######################
function Base.show(io::IO, v::AbstractNormalToricVariety)
    # initiate properties string
    properties_string = ["Normal"]
    
    affine = push_attribute_if_exists!(properties_string, v, :is_affine, "affine")
    
    simplicial_cb!(a, b) = push_attribute_if_exists!(a, b, :is_orbifold, "simplicial")
    push_attribute_if_exists!(properties_string, v, :is_smooth, "smooth"; callback=simplicial_cb!)
    
    projective = nothing
    if isnothing(affine) || !affine
        complete_cb!(a, b) = push_attribute_if_exists!(a, b, :is_complete, "complete")
        projective = push_attribute_if_exists!(properties_string, v, :is_projective, "projective"; callback=complete_cb!)
    end
    
    q_gor_cb!(a, b) = push_attribute_if_exists!(a, b, :is_q_gorenstein, "q-gorenstein")
    gorenstein = push_attribute_if_exists!(properties_string, v, :is_gorenstein, "gorenstein"; callback=q_gor_cb!)
    
    push_attribute_if_exists!(properties_string, v, :is_fano, "fano")
    
    if has_attribute(v, :dim)
        push!(properties_string, string(dim(v))*"-dimensional")
    end
    
    properties_string = [join(properties_string, ", ")]
    push!(properties_string, "toric variety")
    
    push_attribute_if_exists!(properties_string, v, :has_torusfactor, "with torusfactor", "without torusfactor")
    
    join(io, properties_string, " ")
end
