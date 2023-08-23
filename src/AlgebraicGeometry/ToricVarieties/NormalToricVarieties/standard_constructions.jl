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
    variety = normal_toric_variety(C; set_attributes=set_attributes)
    
    if set_attributes
        set_attribute!(variety, :is_complete, false)
        set_attribute!(variety, :is_projective, false)
        set_attribute!(variety, :isprojective_space, false)
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

    pmntv = Polymake.fulton.weighted_projective_space(w)
    variety = NormalToricVariety(pmntv)
    
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
    hirzebruch_surface(::Type{NormalToricVariety}, r::Int; set_attributes::Bool = true)

Constructs the r-th Hirzebruch surface.

# Examples
```jldoctest
julia> hirzebruch_surface(NormalToricVariety, 5)
Normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor
```
"""
function hirzebruch_surface(::Type{NormalToricVariety}, r::Int; set_attributes::Bool = true)
    fan_rays = [1 0; 0 1; -1 r; 0 -1]
    cones = [[1, 2], [2, 3], [3, 4], [4, 1]]
    variety = normal_toric_variety(fan_rays, cones; non_redundant = true)
    
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
    del_pezzo_surface(::Type{NormalToricVariety}, b::Int; set_attributes::Bool = true)

Constructs the del Pezzo surface with `b` blowups for `b` at most 3.

# Examples
```jldoctest
julia> del_pezzo_surface(NormalToricVariety, 3)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor
```
"""
function del_pezzo_surface(::Type{NormalToricVariety}, b::Int; set_attributes::Bool = true)
    # check for valid input
    @req b >= 0 "Number of blowups for construction of del Pezzo surfaces must be non-negative"
    @req b <= 3 "Del Pezzo surfaces with more than three blowups are realized as subvarieties of toric ambient spaces. This is currently not supported"
    
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
    variety = normal_toric_variety(fan_rays, cones; non_redundant = true)
    
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

Blow up the toric variety by subdividing the cone in the list
of *all* cones of the fan of `v` which corresponds to the
provided ideal `I`. Note that this cone need not be maximal.

By default, we pick "e" as the name of the homogeneous coordinate for
the exceptional divisor. As third optional argument one can supply
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

julia> bP3 = domain(blow_up(P3, I))
Normal toric variety

julia> cox_ring(bP3)
Multivariate polynomial ring in 5 variables over QQ graded by
  x1 -> [1 0]
  x2 -> [0 1]
  x3 -> [0 1]
  x4 -> [1 0]
  e -> [1 -1]
```
"""
function blow_up(v::AbstractNormalToricVariety, I::MPolyIdeal; coordinate_name::String = "e", set_attributes::Bool = true)
    @req base_ring(I) == cox_ring(v) "The ideal must be contained in the cox ring of the toric variety"
    indices = [findfirst(y -> y == x, gens(cox_ring(v))) for x in gens(I)]
    @req length(indices) == ngens(I) "All generators must be indeterminates of the cox ring of the toric variety"
    rs = matrix(ZZ, rays(v))
    new_ray = vec(sum([rs[i,:] for i in indices]))
    new_ray = new_ray ./ gcd(new_ray)
    return blow_up(v, new_ray; coordinate_name = coordinate_name, set_attributes = set_attributes)
end


@doc raw"""
    blow_up(v::AbstractNormalToricVariety, new_ray::AbstractVector{<:IntegerUnion}; coordinate_name::String = "e", set_attributes::Bool = true)

Blow up the toric variety by subdividing the fan of the variety with the
provided new ray. Note that this ray must be a primitive element in the
lattice Z^d, with d the dimension of the fan. This function returns the
corresponding blowdown morphism.

By default, we pick "e" as the name of the homogeneous coordinate for
the exceptional divisor. As third optional argument one can supply
a custom variable name.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal, non-affine, smooth, projective, gorenstein, fano, 3-dimensional toric variety without torusfactor

julia> blow_down_morphism = blow_up(P3, [0, 1, 1])
A toric morphism

julia> bP3 = domain(blow_down_morphism)
Normal toric variety

julia> cox_ring(bP3)
Multivariate polynomial ring in 5 variables over QQ graded by
  x1 -> [1 0]
  x2 -> [0 1]
  x3 -> [0 1]
  x4 -> [1 0]
  e -> [1 -1]
```
"""
function blow_up(v::AbstractNormalToricVariety, new_ray::AbstractVector{<:IntegerUnion}; coordinate_name::String = "e", set_attributes::Bool = true)
  return _blow_up(v, star_subdivision(v, new_ray); coordinate_name = coordinate_name, set_attributes = set_attributes)
end



@doc raw"""
    blow_up(v::AbstractNormalToricVariety, n::Int; coordinate_name::String = "e", set_attributes::Bool = true)

Blow up the toric variety by subdividing the n-th cone in the list
of *all* cones of the fan of `v`. This cone need not be maximal.
This function returns the corresponding blowdown morphism.

By default, we pick "e" as the name of the homogeneous coordinate for
the exceptional divisor. As third optional argument one can supply
a custom variable name.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal, non-affine, smooth, projective, gorenstein, fano, 3-dimensional toric variety without torusfactor

julia> blow_down_morphism = blow_up(P3, 5)
A toric morphism

julia> bP3 = domain(blow_down_morphism)
Normal toric variety

julia> cox_ring(bP3)
Multivariate polynomial ring in 5 variables over QQ graded by
  x1 -> [1 0]
  x2 -> [0 1]
  x3 -> [0 1]
  x4 -> [1 0]
  e -> [1 -1]
```
"""
function blow_up(v::AbstractNormalToricVariety, n::Int; coordinate_name::String = "e", set_attributes::Bool = true)
  return _blow_up(v, star_subdivision(v, n); coordinate_name = coordinate_name, set_attributes = set_attributes)
end



function _blow_up(v::AbstractNormalToricVariety, new_fan::PolyhedralFan{QQFieldElem}; coordinate_name::String = "e", set_attributes::Bool = true)
  new_variety = normal_toric_variety(new_fan)
  new_rays = rays(new_fan)
  old_rays = rays(v)
  old_vars = string.(symbols(cox_ring(v)))
  @req !(coordinate_name in old_vars) "The name for the blowup coordinate is already taken"
  new_vars = Vector{String}(undef, length(new_rays))
  for i in 1:length(new_rays)
      j = findfirst(==(new_rays[i]), old_rays)
      new_vars[i] = j !== nothing ? old_vars[j] : coordinate_name
  end
  if set_attributes
    set_attribute!(new_variety, :coordinate_names, new_vars)
  end
  dim = ambient_dim(polyhedral_fan(v))
  return toric_morphism(new_variety, identity_matrix(ZZ, dim), v; check=false)
end



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
Multivariate polynomial ring in 6 variables over QQ graded by
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
Multivariate polynomial ring in 6 variables over QQ graded by
  x1 -> [1 0]
  x2 -> [1 0]
  x3 -> [1 0]
  y1 -> [0 1]
  y2 -> [0 1]
  y3 -> [0 1]
```
"""
function Base.:*(v::AbstractNormalToricVariety, w::AbstractNormalToricVariety; set_attributes::Bool = true)
    product = normal_toric_variety(polyhedral_fan(v)*polyhedral_fan(w))
    set_coordinate_names(product, vcat(["x$(i)" for i in coordinate_names(v)], ["y$(i)" for i in coordinate_names(w)]))
    return product
end



@doc raw"""
    normal_toric_variety_from_star_triangulation(P::Polyhedron; set_attributes::Bool = true)

Returns a toric variety that was obtained from a fine regular
star triangulation of the lattice points of the polyhedron P.
This is particularly useful when the lattice points of the
polyhedron in question admit many triangulations.

# Examples
```jldoctest
julia> P = convex_hull([0 0 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1])
Polyhedron in ambient dimension 3

julia> v = normal_toric_variety_from_star_triangulation(P)
Normal toric variety
```
"""
function normal_toric_variety_from_star_triangulation(P::Polyhedron; set_attributes::Bool = true)
  # Find position of origin in the lattices points of the polyhedron P
  pts = matrix(ZZ, lattice_points(P))
  zero = [0 for i in 1:ambient_dim(P)]
  indices = findall(k -> pts[k,:] == matrix(ZZ, [zero]), 1:nrows(pts))
  @req length(indices) == 1 "Polyhedron must contain origin (exactly once)"

  # Change order of lattice points s.t. zero is the first point
  tmp = pts[1,:]
  pts[1,:] = pts[indices[1],:]
  pts[indices[1],:] = tmp

  # Find one triangulation and turn it into the maximal cones of the toric variety in question. Note that:
  # (a) needs to be converted to incidence matrix
  # (b) one has to remove origin from list of indices (as removed above)
  max_cones = IncidenceMatrix([[c[i]-1 for i in 2:length(c)] for c in _find_full_star_triangulation(pts)])

  # Rays are all but the zero vector at the first position of pts
  integral_rays = vcat([pts[k,:] for k in 2:nrows(pts)])

  # construct the variety
  return normal_toric_variety(integral_rays, max_cones; non_redundant = true)
end



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
ideal(x1*x4)

julia> stanley_reisner_ideal(v2)
ideal(x2*x3)
```
"""
function normal_toric_varieties_from_star_triangulations(P::Polyhedron; set_attributes::Bool = true)
    
    # find position of origin in the lattices points of the polyhedron P
    pts = matrix(ZZ, lattice_points(P))
    zero = [0 for i in 1:ambient_dim(P)]
    indices = findall(k -> pts[k,:] == matrix(ZZ, [zero]), 1:nrows(pts))
    @req length(indices) == 1 "Polyhedron must contain origin (exactly once)"
    
    # change order of lattice points s.t. zero is the first point
    tmp = pts[1,:]
    pts[1,:] = pts[indices[1],:]
    pts[indices[1],:] = tmp
    
    # triangulate the points - note that we enforce full, i.e. want all points to be rays
    trias = star_triangulations(pts; full = true)
    
    # rays are all but the zero vector at the first position of pts
    integral_rays = vcat([pts[k,:] for k in 2:nrows(pts)])
    
    # trias contains the max_cones as list of lists
    # (a) needs to be converted to incidence matrix
    # (b) one has to remove origin from list of indices (as removed above)
    max_cones = [IncidenceMatrix([[c[i]-1 for i in 2:length(c)] for c in t]) for t in trias]
    
    # construct the varieties
    return [normal_toric_variety(integral_rays, cones; non_redundant = true) for cones in max_cones]
end



@doc raw"""
  normal_toric_variety_from_glsm(charges::ZZMatrix; set_attributes::Bool = true)

This function returns one toric variety with the desired
GLSM charges. This can be particularly useful provided that
there are many such toric varieties.

# Examples
```jldoctest
julia> charges = [[1, 1, 1]]
1-element Vector{Vector{Int64}}:
 [1, 1, 1]

julia> normal_toric_variety_from_glsm(charges)
Normal toric variety
```

For convenience, we also support:
- normal_toric_variety_from_glsm(charges::Vector{Vector{Int}})
- normal_toric_variety_from_glsm(charges::Vector{Vector{ZZRingElem}})
"""
function normal_toric_variety_from_glsm(charges::ZZMatrix; set_attributes::Bool = true)

  # find the ray generators
  G1 = free_abelian_group(ncols(charges))
  G2 = free_abelian_group(nrows(charges))
  map = hom(G1, G2, transpose(charges))
  ker = kernel(map)
  embedding = snf(ker[1])[2] * ker[2]
  integral_rays = transpose(embedding.map)

  # identify the points to be triangulated
  pts = matrix(ZZ, zeros(nrows(integral_rays)+1, ncols(integral_rays)))
  pts[2:end, :] = integral_rays

  # construct varieties
  triang = _find_full_star_triangulation(pts)
  max_cones = IncidenceMatrix([[c[i]-1 for i in 2:length(c)] for c in triang])
  variety = normal_toric_variety(integral_rays, max_cones; non_redundant = true)

  # set the attributes and return the variety
  if set_attributes
    set_attribute!(variety, :map_from_torusinvariant_weil_divisor_group_to_class_group, map)
    set_attribute!(variety, :class_group, G2)
    set_attribute!(variety, :torusinvariant_weil_divisor_group, G1)
  end
  return variety
end
normal_toric_variety_from_glsm(charges::Vector{Vector{T}}; set_attributes::Bool = true) where {T <: IntegerUnion} = normal_toric_variety_from_glsm(matrix(ZZ, charges); set_attributes = set_attributes)



@doc raw"""
    normal_toric_varieties_from_glsm(charges::ZZMatrix; set_attributes::Bool = true)

This function returns all toric variety with the desired
GLSM charges. This computation may take a long time if
there are many such toric varieties.

# Examples
```jldoctest
julia> charges = [[1, 1, 1]]
1-element Vector{Vector{Int64}}:
 [1, 1, 1]

julia> normal_toric_varieties_from_glsm(charges)
1-element Vector{NormalToricVariety}:
 Normal toric variety

julia> varieties = normal_toric_varieties_from_glsm(matrix(ZZ, [1 2 3 4 6 0; -1 -1 -2 -2 -3 1]))
1-element Vector{NormalToricVariety}:
 Normal toric variety

julia> cox_ring(varieties[1])
Multivariate polynomial ring in 6 variables over QQ graded by 
  x1 -> [1 -1]
  x2 -> [2 -1]
  x3 -> [3 -2]
  x4 -> [4 -2]
  x5 -> [6 -3]
  x6 -> [0 1]
```

For convenience, we also support:
- normal_toric_varieties_from_glsm(charges::Vector{Vector{Int}})
- normal_toric_varieties_from_glsm(charges::Vector{Vector{ZZRingElem}})
"""
function normal_toric_varieties_from_glsm(charges::ZZMatrix; set_attributes::Bool = true)

    # find the ray generators
    G1 = free_abelian_group(ncols(charges))
    G2 = free_abelian_group(nrows(charges))
    map = hom(G1, G2, transpose(charges))
    ker = kernel(map)
    embedding = snf(ker[1])[2] * ker[2]
    rays = transpose(embedding.map)
    
    # identify the points to be triangulated
    pts = zeros(QQ, nrows(rays), ncols(charges)-nrows(charges))
    for i in 1:nrows(rays)
        pts[i, :] = [ZZRingElem(c) for c in rays[i, :]]
    end
    zero = [0 for i in 1:ncols(charges)-nrows(charges)]
    pts = vcat(matrix(QQ, transpose(zero)), matrix(QQ, pts))
    
    # construct varieties
    integral_rays = vcat([pts[k,:] for k in 2:nrows(pts)])
    max_cones = [IncidenceMatrix([[c[i]-1 for i in 2:length(c)] for c in t]) for t in star_triangulations(pts; full = true)]
    varieties = [normal_toric_variety(integral_rays, cones; non_redundant = true) for cones in max_cones]
    
    # set the map from Div_T -> Cl to the desired matrix
    if set_attributes
      for v in varieties
        set_attribute!(v, :map_from_torusinvariant_weil_divisor_group_to_class_group, map)
        set_attribute!(v, :class_group, G2)
        set_attribute!(v, :torusinvariant_weil_divisor_group, G1)
      end
    end
    
    # return the varieties
    return varieties
end
normal_toric_varieties_from_glsm(charges::Vector{Vector{T}}; set_attributes::Bool = true) where {T <: IntegerUnion} = normal_toric_varieties_from_glsm(matrix(ZZ, charges); set_attributes = set_attributes)
