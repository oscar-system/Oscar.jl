######################
# Julia type for ToricVarieties
######################

pm_object(v::NormalToricVarietyType) = v.polymakeNTV

coefficient_field(::NormalToricVarietyType) = QQ

######################
# Constructors
######################

@doc raw"""
    affine_normal_toric_variety(C::Cone)

Construct the affine normal toric variety $U_{C}$ corresponding to a polyhedral
cone `C`.

# Examples
Set `C` to be the positive orthant in two dimensions.
```jldoctest
julia> C = positive_hull([1 0; 0 1])
Polyhedral cone in ambient dimension 2

julia> antv = affine_normal_toric_variety(C)
Normal toric variety
```
"""
function affine_normal_toric_variety(C::Cone)
  fan = polyhedral_fan(C)
  pmntv = Polymake.fulton.NormalToricVariety(Oscar.pm_object(fan))
  variety = AffineNormalToricVariety(pmntv)
  set_attribute!(variety, :cone, C)
  set_attribute!(variety, :fan, fan)
  return variety
end


@doc raw"""
    normal_toric_variety(C::Cone)

Construct the (affine) normal toric variety $X_{\Sigma}$ corresponding to a
polyhedral fan $\Sigma = C$ consisting only of the cone `C`.

# Examples
Set `C` to be the positive orthant in two dimensions.
```jldoctest
julia> C = positive_hull([1 0; 0 1])
Polyhedral cone in ambient dimension 2

julia> ntv = normal_toric_variety(C)
Normal toric variety
```
"""
function normal_toric_variety(C::Cone)
  fan = polyhedral_fan(C)
  pmntv = Polymake.fulton.NormalToricVariety(Oscar.pm_object(fan))
  variety = NormalToricVariety(pmntv)
  set_attribute!(variety, :cone, C)
  set_attribute!(variety, :fan, fan)
  return variety
end


@doc raw"""
    normal_toric_variety(max_cones::IncidenceMatrix, rays::AbstractCollection[RayVector]; non_redundant::Bool = false)

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

julia> max_cones = IncidenceMatrix([[1, 2], [2, 3], [3, 4], [4, 1]])
4×4 IncidenceMatrix
[1, 2]
[2, 3]
[3, 4]
[1, 4]

julia> normal_toric_variety(max_cones, ray_generators)
Normal toric variety

julia> normal_toric_variety(max_cones, ray_generators; non_redundant = true)
Normal toric variety
```
"""
normal_toric_variety(max_cones::Vector{Vector{Int64}}, rays::AbstractCollection[RayVector]; non_redundant::Bool = false) = normal_toric_variety(IncidenceMatrix(max_cones), rays; non_redundant)
function normal_toric_variety(max_cones::IncidenceMatrix, rays::AbstractCollection[RayVector]; non_redundant::Bool = false)
  fan = polyhedral_fan(max_cones, rays; non_redundant)
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
  return NormalToricVariety(pmntv)
end


@doc raw"""
    normal_toric_variety(P::Polyhedron)

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
function normal_toric_variety(P::Polyhedron)
  variety = normal_toric_variety(normal_fan(P))
  set_attribute!(variety, :polyhedron, P)
  return variety
end


@doc raw"""
    affine_normal_toric_variety(v::NormalToricVariety)

For internal design, we make a strict distinction between
normal toric varieties and affine toric varieties.
Given an affine, normal toric variety `v`,
this method turns it into an affine toric variety.

# Examples
```jldoctest
julia> v = normal_toric_variety(positive_hull([1 0; 0 1]))
Normal toric variety

julia> affineVariety = affine_normal_toric_variety(v)
Normal toric variety
```
"""
function affine_normal_toric_variety(v::NormalToricVariety)
  is_affine(v) || error("Cannot construct affine toric variety from non-affine input")
  return AffineNormalToricVariety(pm_object(v))
end



######################
# Display
######################

function Base.show(io::IO, v::NormalToricVarietyType)
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

Base.show(io::IO, ::MIME"text/plain", v::NormalToricVarietyType) = Base.show(pretty(io), v)
