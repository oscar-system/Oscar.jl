######################
# 1: The Julia type for ToricVarieties
######################
abstract type AbstractNormalToricVariety <: _FanLikeType{QQFieldElem} end

@attributes mutable struct NormalToricVariety <: AbstractNormalToricVariety
           polymakeNTV::Polymake.BigObject
           NormalToricVariety(polymakeNTV::Polymake.BigObject) = new(polymakeNTV)
end

@attributes mutable struct AffineNormalToricVariety <: AbstractNormalToricVariety
           polymakeNTV::Polymake.BigObject
           AffineNormalToricVariety(polymakeNTV::Polymake.BigObject) = new(polymakeNTV)
end

pm_object(v::AbstractNormalToricVariety) = v.polymakeNTV

######################
# 2: Generic constructors
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
    fan = PolyhedralFan(C)
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
    fan = PolyhedralFan(C)
    pmntv = Polymake.fulton.NormalToricVariety(Oscar.pm_object(fan))
    variety = NormalToricVariety(pmntv)
    
    if set_attributes
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
    normal_toric_variety(rays::Vector{Vector{Int64}}, max_cones::Vector{Vector{Int64}}; non_redundant::Bool = false, set_attributes::Bool = true)

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
function normal_toric_variety(rays::Vector{Vector{Int64}}, max_cones::Vector{Vector{Int64}}; non_redundant::Bool = false, set_attributes::Bool = true)
    fan = PolyhedralFan(transpose(hcat(rays...)), IncidenceMatrix(max_cones); non_redundant = non_redundant)
    return normal_toric_variety(fan; set_attributes = set_attributes)
end


@doc raw"""
    normal_toric_variety(PF::PolyhedralFan; set_attributes::Bool = true)

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
function normal_toric_variety(PF::PolyhedralFan; set_attributes::Bool = true)
    fan = Oscar.pm_object(PF)
    pmntv = Polymake.fulton.NormalToricVariety(fan)
    variety = NormalToricVariety(pmntv)
    if set_attributes
        set_attribute!(variety, :fan, PF)
    end
    return variety
end


@doc raw"""
    NormalToricVariety(P::Polyhedron; set_attributes::Bool = true)

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
        set_attribute!(variety, :is_complete, false)
        set_attribute!(variety, :is_projective, false)
        set_attribute!(variety, :is_projective_space, false)
    end
    
    return variety
end


######################
# 3: Special constructors
######################

@doc raw"""
    affine_space(::Type{NormalToricVariety}, d::Int; set_attributes::Bool = true)

Constructs the (toric) affine space of dimension `d`.

# Examples
```jldoctest
julia> affine_space(NormalToricVariety, 2)
Normal, affine, 2-dimensional toric variety
```
"""
function affine_space(::Type{NormalToricVariety}, d::Int; set_attributes::Bool = true)
    C = positive_hull(identity_matrix(ZZ, d))
    fan = PolyhedralFan(C)
    pmntv = Polymake.fulton.NormalToricVariety(Oscar.pm_object(fan))
    variety = NormalToricVariety(pmntv)
    
    if set_attributes
        set_attribute!(variety, :is_affine, true)
        set_attribute!(variety, :is_complete, false)
        set_attribute!(variety, :is_projective, false)
        set_attribute!(variety, :isprojective_space, false)
        set_attribute!(variety, :fan, fan)
        set_attribute!(variety, :dim, d)
        set_attribute!(variety, :dim_of_torusfactor, 0)
    end
    
    return variety
end


@doc raw"""
    projective_space(::Type{NormalToricVariety}, d::Int; set_attributes::Bool = true)

Construct the projective space of dimension `d`.

# Examples
```jldoctest
julia> projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor
```
"""
function projective_space(::Type{NormalToricVariety}, d::Int; set_attributes::Bool = true)
    f = normal_fan(Oscar.simplex(d))
    pm_object = Polymake.fulton.NormalToricVariety(Oscar.pm_object(f))
    variety = NormalToricVariety(pm_object)
    
    # make standard choice for weights of the cox ring
    set_attribute!(variety, :torusinvariant_cartier_divisor_group, free_abelian_group(d+1))
    set_attribute!(variety, :torusinvariant_weil_divisor_group, free_abelian_group(d+1))
    set_attribute!(variety, :class_group, free_abelian_group(1))
    set_attribute!(variety, :picard_group, free_abelian_group(1))
    weights = matrix(ZZ,ones(Int,d+1,1))
    set_attribute!(variety, :map_from_torusinvariant_weil_divisor_group_to_class_group, hom(torusinvariant_weil_divisor_group(variety), class_group(variety), weights))
    set_attribute!(variety, :map_from_torusinvariant_cartier_divisor_group_to_picard_group, hom(torusinvariant_cartier_divisor_group(variety), picard_group(variety), weights))
    set_attribute!(variety, :map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group, hom(torusinvariant_cartier_divisor_group(variety), torusinvariant_weil_divisor_group(variety), identity_matrix(ZZ,d+1)))
    
    if set_attributes
        set_attribute!(variety, :is_affine, false)
        set_attribute!(variety, :is_projective, true)
        set_attribute!(variety, :is_projective_space, true)
        set_attribute!(variety, :is_smooth, true)
        set_attribute!(variety, :is_complete, true)
        set_attribute!(variety, :has_torusfactor, false)
        set_attribute!(variety, :is_orbifold, true)
        set_attribute!(variety, :is_simplicial, true)
        set_attribute!(variety, :is_gorenstein, true)
        set_attribute!(variety, :is_q_gorenstein, true)
        set_attribute!(variety, :is_fano, true)
        set_attribute!(variety, :dim, d)
        set_attribute!(variety, :dim_of_torusfactor, 0)
        set_attribute!(variety, :euler_characteristic, d+1)
        set_attribute!(variety, :betti_number, fill(ZZRingElem(1), d+1))
        set_attribute!(variety, :character_lattice, free_abelian_group(d))
    end
    
    return variety
end


@doc raw"""
    weighted_projective_space(::Type{NormalToricVariety}, w::Vector{T}; set_attributes::Bool = true) where {T <: IntegerUnion}

Construct the weighted projective space corresponding to the weights `w`.

# Examples
```jldoctest
julia> weighted_projective_space(NormalToricVariety, [2,3,1])
Normal, non-affine, simplicial, projective, 2-dimensional toric variety without torusfactor
```
"""
function weighted_projective_space(::Type{NormalToricVariety}, w::Vector{T}; set_attributes::Bool = true) where {T <: IntegerUnion}
    # build projective space
    if all(a -> isone(a), w)
      return projective_space(NormalToricVariety, length(w)-1)
    end

    # construct the weighted projective space
    # pmntv = Polymake.fulton.weighted_projective_space(w)
    # variety = NormalToricVariety(pmntv)
    
    # Workaround due to bug in polymake
    #
    # We follow the recipe of page 35-36 in
    # William Fulton: Introduction to toric varieties
    # Since we cannot deal with lattices that are finer than ZZ^n, we scale
    # everything up by the corresponding lcms.
    lcms = [lcm(w[1], w[i]) for i in 2:length(w)]
    lattice_gens = ZZMatrix(length(w), length(w)-1)
    ray_gens = ZZMatrix(length(w), length(w)-1)
    lattice_gens[1,:] = [-div(i, w[1]) for i in lcms]
    ray_gens[1,:] = -lcms
    for i in 1:length(w)-1
        lattice_gens[i+1,i] = div(lcms[i], w[i+1])
        ray_gens[i+1,i] = lcms[i]
    end
    H = hnf(lattice_gens)
    lattice_gens = transpose(H[1:length(w)-1, :])
    tr,_ = pseudo_inv(lattice_gens)
    ray_gens = ray_gens * transpose(tr)
    mc = IncidenceMatrix(subsets(Vector{Int}(1:length(w)), length(w)-1))
    variety = normal_toric_variety(PolyhedralFan(ray_gens, mc; non_redundant=true ))
    
    # make standard choice for the weights of the cox ring
    set_attribute!(variety, :torusinvariant_weil_divisor_group, free_abelian_group(length(w)))
    set_attribute!(variety, :class_group, free_abelian_group(1))
    weights = matrix(ZZ,hcat(w))
    set_attribute!(variety, :map_from_torusinvariant_weil_divisor_group_to_class_group, hom(torusinvariant_weil_divisor_group(variety), class_group(variety), weights))
    
    if set_attributes
        set_attribute!(variety, :has_torusfactor, false)
        set_attribute!(variety, :is_affine, false)
        set_attribute!(variety, :is_projective, true)
        set_attribute!(variety, :is_projective_space, false)
        set_attribute!(variety, :is_complete, true)
        set_attribute!(variety, :is_orbifold, true)
        set_attribute!(variety, :is_simplicial, true)
        set_attribute!(variety, :dim, length(w)-1)
        set_attribute!(variety, :character_lattice, free_abelian_group(length(w)-1))
        set_attribute!(variety, :dim_of_torusfactor, 0)
    end
    
    return variety
end


@doc raw"""
    hirzebruch_surface(r::Int; set_attributes::Bool = true)

Constructs the r-th Hirzebruch surface.

# Examples
```jldoctest
julia> hirzebruch_surface(5)
Normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor
```
"""
function hirzebruch_surface(r::Int; set_attributes::Bool = true)
    fan_rays = [1 0; 0 1; -1 r; 0 -1]
    cones = IncidenceMatrix([[1, 2], [2, 3], [3, 4], [4, 1]])
    variety = normal_toric_variety(PolyhedralFan(fan_rays, cones; non_redundant = true))
    
    # make standard choice for the weights of the cox ring
    set_attribute!(variety, :torusinvariant_cartier_divisor_group, free_abelian_group(4))
    set_attribute!(variety, :torusinvariant_weil_divisor_group, free_abelian_group(4))
    set_attribute!(variety, :class_group, free_abelian_group(2))
    set_attribute!(variety, :picard_group, free_abelian_group(2))
    weights = matrix(ZZ, [1 0; 0 1; 1 0; r 1])
    set_attribute!(variety, :map_from_torusinvariant_weil_divisor_group_to_class_group, hom(torusinvariant_weil_divisor_group(variety), class_group(variety), weights))
    set_attribute!(variety, :map_from_torusinvariant_cartier_divisor_group_to_picard_group, hom(torusinvariant_cartier_divisor_group(variety), picard_group(variety), weights))
    set_attribute!(variety, :map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group, hom(torusinvariant_cartier_divisor_group(variety), torusinvariant_weil_divisor_group(variety), identity_matrix(ZZ,4)))
    
    if set_attributes
        set_attribute!(variety, :is_affine, false)
        set_attribute!(variety, :is_projective, true)
        set_attribute!(variety, :is_projective_space, false)
        set_attribute!(variety, :is_smooth, true)
        set_attribute!(variety, :is_complete, true)
        set_attribute!(variety, :has_torusfactor, false)
        set_attribute!(variety, :is_orbifold, true)
        set_attribute!(variety, :is_simplicial, true)
        set_attribute!(variety, :is_gorenstein, true)
        set_attribute!(variety, :is_q_gorenstein, true)
        if abs(r) <= 1
            set_attribute!(variety, :is_fano, true)
        else
            set_attribute!(variety, :is_fano, false)
        end
        vars = ["t1", "x1", "t2", "x2"]
        set_coordinate_names(variety, vars)
        set_attribute!(variety, :dim, 2)
        set_attribute!(variety, :dim_of_torusfactor, 0)
        set_attribute!(variety, :euler_characteristic, 4)
        set_attribute!(variety, :betti_number, [ZZRingElem(1), ZZRingElem(2), ZZRingElem(1)])
        set_attribute!(variety, :character_lattice, free_abelian_group(2))
        set_attribute!(variety, :map_from_character_lattice_to_torusinvariant_weil_divisor_group, hom(character_lattice(variety), torusinvariant_weil_divisor_group(variety), transpose(matrix(ZZ, fan_rays))))
    end
    
    return variety
end


@doc raw"""
    del_pezzo_surface(b::Int; set_attributes::Bool = true)

Constructs the del Pezzo surface with `b` blowups for `b` at most 3.

# Examples
```jldoctest
julia> del_pezzo_surface(3)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor
```
"""
function del_pezzo_surface(b::Int; set_attributes::Bool = true)
    # check for valid input
    if b < 0
        throw(ArgumentError("Number of blowups for construction of del Pezzo surfaces must be non-negative"))
    end
    if b > 3
        throw(ArgumentError("Del Pezzo surfaces with more than three blowups are realized as subvarieties of toric ambient spaces. This is currently not supported"))
    end
    
    # special case of projective space
    if b == 0
        return projective_space(NormalToricVariety, 2)
    end
    
    # construct the "true" toric del Pezzo surfaces
    if b == 1
        fan_rays = [1 0; 0 1; -1 -1; 1 1]
        cones = IncidenceMatrix([[1, 4], [2, 4], [1, 3], [2, 3]])
    end
    if b == 2
        fan_rays = [1 0; 0 1; -1 -1; 1 1; 0 -1]
        cones = IncidenceMatrix([[1, 4], [2, 4], [1, 5], [5, 3], [2, 3]])
    end
    if b == 3
        fan_rays = [1 0; 0 1; -1 -1; 1 1; 0 -1; -1 0]
        cones = IncidenceMatrix([[1, 4], [2, 4], [1, 5], [5, 3], [2, 6], [6, 3]])
    end
    variety = normal_toric_variety(PolyhedralFan(fan_rays, cones; non_redundant = true))
    
    # make standard choice for weights of the cox ring
    if b == 1
        set_attribute!(variety, :torusinvariant_cartier_divisor_group, free_abelian_group(4))
        set_attribute!(variety, :torusinvariant_weil_divisor_group, free_abelian_group(4))
        set_attribute!(variety, :class_group, free_abelian_group(2))
        set_attribute!(variety, :picard_group, free_abelian_group(2))
        weights = matrix(ZZ, [1 1; 1 1; 1 0; 0 -1])
        set_attribute!(variety, :map_from_torusinvariant_weil_divisor_group_to_class_group, hom(torusinvariant_weil_divisor_group(variety), class_group(variety), weights))
        set_attribute!(variety, :map_from_torusinvariant_cartier_divisor_group_to_picard_group, hom(torusinvariant_cartier_divisor_group(variety), picard_group(variety), weights))
        set_attribute!(variety, :map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group, hom(torusinvariant_cartier_divisor_group(variety), torusinvariant_weil_divisor_group(variety), identity_matrix(ZZ,4)))
    end
    if b == 2
        set_attribute!(variety, :torusinvariant_cartier_divisor_group, free_abelian_group(5))
        set_attribute!(variety, :torusinvariant_weil_divisor_group, free_abelian_group(5))
        set_attribute!(variety, :class_group, free_abelian_group(3))
        set_attribute!(variety, :picard_group, free_abelian_group(3))
        weights = matrix(ZZ, [1 1 1; 1 1 0; 1 0 1; 0 -1 0; 0 0 -1])
        set_attribute!(variety, :map_from_torusinvariant_weil_divisor_group_to_class_group, hom(torusinvariant_weil_divisor_group(variety), class_group(variety), weights))
        set_attribute!(variety, :map_from_torusinvariant_cartier_divisor_group_to_picard_group, hom(torusinvariant_cartier_divisor_group(variety), picard_group(variety), weights))
        set_attribute!(variety, :map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group, hom(torusinvariant_cartier_divisor_group(variety), torusinvariant_weil_divisor_group(variety), identity_matrix(ZZ,5)))
    end
    if b == 3
        set_attribute!(variety, :torusinvariant_cartier_divisor_group, free_abelian_group(6))
        set_attribute!(variety, :torusinvariant_weil_divisor_group, free_abelian_group(6))
        set_attribute!(variety, :class_group, free_abelian_group(4))
        set_attribute!(variety, :picard_group, free_abelian_group(4))
        weights = matrix(ZZ, [1 1 1 0; 1 1 0 1; 1 0 1 1; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1])
        set_attribute!(variety, :map_from_torusinvariant_weil_divisor_group_to_class_group, hom(torusinvariant_weil_divisor_group(variety), class_group(variety), weights))
        set_attribute!(variety, :map_from_torusinvariant_cartier_divisor_group_to_picard_group, hom(torusinvariant_cartier_divisor_group(variety), picard_group(variety), weights))
        set_attribute!(variety, :map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group, hom(torusinvariant_cartier_divisor_group(variety), torusinvariant_weil_divisor_group(variety), identity_matrix(ZZ,6)))
    end
    
    if set_attributes
        set_attribute!(variety, :is_affine, false)
        set_attribute!(variety, :is_projective, true)
        set_attribute!(variety, :is_projective_space, false)
        set_attribute!(variety, :is_smooth, true)
        set_attribute!(variety, :is_complete, true)
        set_attribute!(variety, :has_torusfactor, false)
        set_attribute!(variety, :is_orbifold, true)
        set_attribute!(variety, :is_simplicial, true)
        set_attribute!(variety, :is_gorenstein, true)
        set_attribute!(variety, :is_q_gorenstein, true)
        set_attribute!(variety, :is_fano, true)
        vars = ["x1", "x2", "x3", "e1", "e2", "e3"]
        set_coordinate_names(variety, vars[1:(3 + b)])
        set_attribute!(variety, :dim, 2)
        set_attribute!(variety, :dim_of_torusfactor, 0)
        if b == 1
            set_attribute!(variety, :euler_characteristic, 4)
            set_attribute!(variety, :betti_number, [ZZRingElem(1), ZZRingElem(2), ZZRingElem(1)])
            set_attribute!(variety, :character_lattice, free_abelian_group(2))
            set_attribute!(variety, :map_from_character_lattice_to_torusinvariant_weil_divisor_group, hom(character_lattice(variety), torusinvariant_weil_divisor_group(variety), transpose(matrix(ZZ, fan_rays))))
        end
        if b == 2
            set_attribute!(variety, :euler_characteristic, 5)
            set_attribute!(variety, :betti_number, [ZZRingElem(1), ZZRingElem(3), ZZRingElem(1)])
            set_attribute!(variety, :character_lattice, free_abelian_group(2))
            set_attribute!(variety, :map_from_character_lattice_to_torusinvariant_weil_divisor_group, hom(character_lattice(variety), torusinvariant_weil_divisor_group(variety), transpose(matrix(ZZ, fan_rays))))
        end
        if b == 3
            set_attribute!(variety, :euler_characteristic, 6)
            set_attribute!(variety, :betti_number, [ZZRingElem(1), ZZRingElem(4), ZZRingElem(1)])
            set_attribute!(variety, :character_lattice, free_abelian_group(2))
            set_attribute!(variety, :map_from_character_lattice_to_torusinvariant_weil_divisor_group, hom(character_lattice(variety), torusinvariant_weil_divisor_group(variety), transpose(matrix(ZZ, fan_rays))))
        end
    end
    
    return variety
end


############################
# 4: Blowups
############################

@doc raw"""
  blow_up(v::AbstractNormalToricVariety, I::MPolyIdeal; coordinate_name::String = "e", set_attributes::Bool = true)

Blowup the toric variety by subdividing the cone in the list
of *all* cones of the fan of `v` which corresponds to the
provided ideal `I`. Note that this cone need not be maximal.

By default, we pick "e" as the name of the homogeneous coordinate for
the exceptional divisor. As third optional argument on can supply
a custom variable name.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal, non-affine, smooth, projective, gorenstein, fano, 3-dimensional toric variety without torusfactor

julia> (x1,x2,x3,x4) = gens(cox_ring(P3))
4-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3
 x4

julia> I = ideal([x2,x3])
ideal(x2, x3)

julia> bP3 = blow_up(P3, I)
Normal toric variety

julia> cox_ring(bP3)
Multivariate Polynomial Ring in x1, x2, e, x3, x4 over Rational Field graded by 
  x1 -> [1 0]
  x2 -> [0 1]
  e -> [1 -1]
  x3 -> [0 1]
  x4 -> [1 0]
```
"""
function blow_up(v::AbstractNormalToricVariety, I::MPolyIdeal; coordinate_name::String = "e", set_attributes::Bool = true)
    @req parent(gens(I)[1]) == cox_ring(v) "The ideal must be contained in the cox ring of the toric variety"
    indices = [findfirst(y -> y == x, gens(cox_ring(v))) for x in gens(I)]
    @req length(indices) == ngens(I) "All generators must be indeterminates of the cox ring of the toric variety"
    cone_list = cones(v)
    cone_indices = [filter(l -> cone_list[k,l], 1:ncols(cone_list)) for k in 1:n_cones(v)]
    cone_index = findfirst(x -> x == indices, cone_indices)
    @req cone_index !== nothing "There is no corresponding cone that could be subdivided"
    return blow_up(v, cone_index; coordinate_name = coordinate_name, set_attributes = set_attributes)
end


@doc raw"""
  blow_up(v::AbstractNormalToricVariety, n::Int; coordinate_name::String = "e", set_attributes::Bool = true)

Blowup the toric variety by subdividing the n-th cone in the list
of *all* cones of the fan of `v`. This cone need not be maximal.

By default, we pick "e" as the name of the homogeneous coordinate for
the exceptional divisor. As third optional argument on can supply
a custom variable name.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal, non-affine, smooth, projective, gorenstein, fano, 3-dimensional toric variety without torusfactor

julia> cones(P3)
14×4 IncidenceMatrix
[1, 2, 3]
[2, 3, 4]
[1, 3, 4]
[1, 2, 4]
[2, 3]
[1, 3]
[1, 2]
[3, 4]
[2, 4]
[1, 4]
[2]
[3]
[1]
[4]

julia> bP3 = blow_up(P3, 5)
Normal toric variety

julia> cox_ring(bP3)
Multivariate Polynomial Ring in x1, x2, e, x3, x4 over Rational Field graded by 
  x1 -> [1 0]
  x2 -> [0 1]
  e -> [1 -1]
  x3 -> [0 1]
  x4 -> [1 0]
```
"""
function blow_up(v::AbstractNormalToricVariety, n::Int; coordinate_name::String = "e", set_attributes::Bool = true)
    new_fan = starsubdivision(fan(v), n)
    new_variety = normal_toric_variety(new_fan; set_attributes = set_attributes)
    new_rays = rays(new_fan)
    old_vars = [string(x) for x in gens(cox_ring(v))]
    @req findfirst(x->occursin(coordinate_name, x), old_vars) === nothing "The name for the blowup coordinate is already taken"
    new_vars = [if new_rays[i] in rays(fan(v)) old_vars[findfirst(x->x==new_rays[i], rays(fan(v)))] else coordinate_name end for i in 1:length(new_rays)]
    set_attribute!(new_variety, :coordinate_names, new_vars)
    return new_variety
end



############################
# 5: Cartesian product
############################

@doc raw"""
    Base.:*(v::AbstractNormalToricVariety, w::AbstractNormalToricVariety; set_attributes::Bool = true)

Return the Cartesian/direct product of two normal toric varieties `v` and `w`.

By default, we prepend an "x" to all homogeneous coordinate names of the first factor
`v` and a "y" to all homogeneous coordinate names of the second factor `w`. This default
can be overwritten by invoking `set_coordinate_names` after creating the variety
(cf. [`set_coordinate_names(v::AbstractNormalToricVariety, coordinate_names::Vector{String})`](@ref)).

*Important*: Recall that the coordinate names can only be changed as long as the toric
variety in question is not finalized (cf. [`is_finalized(v::AbstractNormalToricVariety)`](@ref)).

Crucially, the order of the homogeneous coordinates is not shuffled. To be more
specific, assume that `v` has ``n_1`` and `w` has ``n_2`` homogeneous coordinates. Then
`v * w` has ``n_1 + n_2`` homogeneous coordinates. The first ``n_1`` of these coordinates
are those of `v` and appear in the very same order as they do for `v`. The remaining
``n_2`` homogeneous coordinates are those of `w` and appear in the very same order as they
do for `w`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> v1 = P2 * P2
Normal toric variety

julia> cox_ring(v1)
Multivariate Polynomial Ring in 6 variables xx1, xx2, xx3, yx1, ..., yx3 over Rational Field graded by
  xx1 -> [1 0]
  xx2 -> [1 0]
  xx3 -> [1 0]
  yx1 -> [0 1]
  yx2 -> [0 1]
  yx3 -> [0 1]

julia> v2 = P2 * P2
Normal toric variety

julia> set_coordinate_names(v2, ["x1", "x2", "x3", "y1", "y2", "y3"])

julia> cox_ring(v2)
Multivariate Polynomial Ring in 6 variables x1, x2, x3, y1, ..., y3 over Rational Field graded by
  x1 -> [1 0]
  x2 -> [1 0]
  x3 -> [1 0]
  y1 -> [0 1]
  y2 -> [0 1]
  y3 -> [0 1]
```
"""
function Base.:*(v::AbstractNormalToricVariety, w::AbstractNormalToricVariety; set_attributes::Bool = true)
    product = normal_toric_variety(fan(v)*fan(w); set_attributes = set_attributes)
    set_coordinate_names(product, vcat(["x$(i)" for i in coordinate_names(v)], ["y$(i)" for i in coordinate_names(w)]))
    return product
end



########################################
# 6: Toric varieties from triangulations
########################################

@doc raw"""
    normal_toric_varieties_from_star_triangulations(P::Polyhedron; set_attributes::Bool = true)

Returns the list of toric varieties obtained from fine regular
star triangulations of the polyhedron P. With this we can
compute the two phases of the famous conifold transition.

# Examples
```jldoctest
julia> P = convex_hull([0 0 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1])
Polyhedron in ambient dimension 3

julia> (v1, v2) = normal_toric_varieties_from_star_triangulations(P)
2-element Vector{NormalToricVariety}:
 Normal toric variety
 Normal toric variety

julia> stanley_reisner_ideal(v1)
ideal(x2*x4)

julia> stanley_reisner_ideal(v2)
ideal(x1*x3)
```
"""
function normal_toric_varieties_from_star_triangulations(P::Polyhedron; set_attributes::Bool = true)
    # triangulate the polyhedron
    trias = star_triangulations(P)
    
    # Currently, the rays in trias[1]
    # (a) are encoded as QQMatrix (ZZRingElem expected)
    # (b) contain the origin as first element (not a rays, so to be removed)
    rays = trias[1]
    integral_rays = zeros(ZZ, nrows(rays)-1, ncols(rays))
    for i in 2:nrows(rays)
        integral_rays[i-1, 1:ncols(rays)] = [ZZ(c) for c in rays[i, 1:ncols(rays)]]
    end
    
    # trias[2] contains the max_cones as list of lists
    # (a) needs to be converted to incidence matrix
    # (b) one has to remove origin from list of indices (as removed above)
    max_cones = trias[2]
    max_cones = [IncidenceMatrix([[c[i]-1 for i in 2:length(c)] for c in t]) for t in max_cones]
    
    # construct the varieties
    return [normal_toric_variety(PolyhedralFan(integral_rays, cones; non_redundant = true), set_attributes = set_attributes) for cones in max_cones]
end



###############################
# 7: Toric varieties from GLSMs
###############################

@doc raw"""
    normal_toric_varieties_from_glsm(charges::ZZMatrix; set_attributes::Bool = true)

Witten's Generalized-Sigma models (GLSM) [Wit88](@cite)
originally sparked interest in the physics community in toric varieties.
On a mathematical level, this establishes a construction of toric
varieties for  which a Z^n grading of the Cox ring is provided. See
for example [FJR17](@cite), which describes this as GIT
construction [CLS11](@cite).

Explicitly, given the grading of the Cox ring, the map from
the group of torus invariant Weil divisors to the class group
is known. Under the assumption that the variety in question
has no torus factor, we can then identify the map from the
lattice to the group of torus invariant Weil divisors as the
kernel of the map from the torus invariant Weil divisor to the
class group. The latter is a map between free Abelian groups, i.e.
is provided by an integer valued matrix. The rows of this matrix
are nothing but the ray generators of the fan of the toric variety.
It then remains to triangulate these rays, hence in general for
a GLSM the toric variety is only unique up to fine regular
star triangulations.

# Examples
```jldoctest
julia> charges = [[1, 1, 1]]
1-element Vector{Vector{Int64}}:
 [1, 1, 1]

julia> normal_toric_varieties_from_glsm(charges)
1-element Vector{NormalToricVariety}:
 Normal toric variety
```

For convenience, we also support:
- normal_toric_varieties_from_glsm(charges::Vector{Vector{Int}})
- normal_toric_varieties_from_glsm(charges::Vector{Vector{ZZRingElem}})
"""
function normal_toric_varieties_from_glsm(charges::ZZMatrix; set_attributes::Bool = true)
    # compute the map from Div_T -> Cl
    source = free_abelian_group(ncols(charges))
    range = free_abelian_group(nrows(charges))
    map = hom(source, range, transpose(charges))
    
    # compute the map char -> Div_T
    ker = kernel(map)
    embedding = snf(ker[1])[2] * ker[2]
    
    # identify the rays
    rays = transpose(embedding.map)
    
    # construct vertices of polyhedron
    pts = zeros(QQ, nrows(rays), ncols(charges)-nrows(charges))
    for i in 1:nrows(rays)
        pts[i, :] = [ZZRingElem(c) for c in rays[i, :]]
    end
    zero = [0 for i in 1:ncols(charges)-nrows(charges)]
    pts = vcat(matrix(QQ, transpose(zero)), matrix(QQ, pts))
    
    # construct polyhedron
    p = convex_hull(pts)
    return normal_toric_varieties_from_star_triangulations(p; set_attributes = set_attributes)
end
normal_toric_varieties_from_glsm(charges::Vector{Vector{T}}; set_attributes::Bool = true) where {T <: IntegerUnion} = normal_toric_varieties_from_glsm(matrix(ZZ, charges); set_attributes = set_attributes)



############################
### 8: Display
############################
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
