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
as vector of vectors. By default, this method allows redundancies in the input, e.g. duplicate rays and non-maximal cones. If the user
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

julia> max_cones = incidence_matrix([[1, 2], [2, 3], [3, 4], [4, 1]])
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
normal_toric_variety(
  max_cones::Vector{Vector{Int64}},
  rays::AbstractCollection[RayVector];
  non_redundant::Bool=false,
) = normal_toric_variety(
  IncidenceMatrix(max_cones), rays; non_redundant
)
function normal_toric_variety(
  max_cones::IncidenceMatrix, rays::AbstractCollection[RayVector]; non_redundant::Bool=false
)
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
Polytope in ambient dimension 2

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
Polytope in ambient dimension 2

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
# Equality
######################

@doc raw"""
    (==)(X::NormalToricVariety, Y::NormalToricVariety) -> Bool

Check equality of the polyhedral fans as sets of cones.

# Examples
```jldoctest
julia> H = hirzebruch_surface(NormalToricVariety, 0)
Normal toric variety

julia> P1 = projective_space(NormalToricVariety, 1)
Normal toric variety

julia> H == P1 * P1
true
```
"""
function Base.:(==)(X::NormalToricVariety, Y::NormalToricVariety)
  X === Y && return true
  ambient_dim(X) == ambient_dim(Y) || return false
  n_rays(X) == n_rays(Y) || return false

  # p is a permutation such that the i-th ray of X is the p(i)-th ray of Y
  p = inv(perm(sortperm(rays(X)))) * perm(sortperm(rays(Y)))

  for i in 1:n_rays(X)
    rays(X)[i] == rays(Y)[p(i)] || return false
  end
  @inline rows(Z) = [
    row(maximal_cones(IncidenceMatrix, Z), i) for i in 1:n_maximal_cones(Z)
  ]
  return Set(map(r -> Set(p.(r)), rows(X))) == Set(rows(Y))
end

@doc raw"""
    _id(X::NormalToricVariety)
    -> Tuple{Vector{Vector{QQFieldElem}}, Vector{Vector{Int64}}}

Given a toric variety `X`, return a pair `Oscar._id(X)` with the
following property: two toric varieties `X` and `Y` have equal
polyhedral fans, taken as sets of cones, if and only if
`Oscar._id(X) == Oscar._id(Y)`.

# Examples
```jldoctest
julia> H = hirzebruch_surface(NormalToricVariety, 0)
Normal toric variety

julia> P1 = projective_space(NormalToricVariety, 1)
Normal toric variety

julia> Oscar._id(H) == Oscar._id(P1 * P1)
true
```
"""
function _id(X::NormalToricVariety)
  p = inv(perm(sortperm(rays(X))))
  sorted_rays = Vector.(permuted(collect(rays(X)), p))
  @inline rows(Z) = [
    row(maximal_cones(IncidenceMatrix, Z), i) for i in 1:n_maximal_cones(Z)
  ]
  sorted_maximal_cones = sort(map(r -> sort(Vector(p.(r))), rows(X)))
  return (sorted_rays, sorted_maximal_cones)
end

function Base.hash(X::NormalToricVariety, h::UInt)
  return hash(_id(X), h)
end

######################
# Display
######################

function Base.show(io::IO, v::NormalToricVarietyType)
  # initiate properties string
  properties_string = []

  if has_attribute(v, :is_smooth) && get_attribute(v, :is_smooth)
    push!(properties_string, "Smooth")
  else
    if has_attribute(v, :is_orbifold) && get_attribute(v, :is_orbifold)
      push!(properties_string, "Q-factorial")
    end

    if has_attribute(v, :is_gorenstein) && get_attribute(v, :is_gorenstein)
      push!(properties_string, "Gorenstein")
    elseif (
      has_attribute(v, :is_q_gorenstein) && get_attribute(v, :is_q_gorenstein)
      && (!has_attribute(v, :is_orbifold) || !get_attribute(v, :is_orbifold))
    )
      push!(properties_string, "Q-Gorenstein")
    end

    if isempty(properties_string)
      push!(properties_string, "Normal")
    else
      push!(properties_string, "normal")
    end
  end

  if has_attribute(v, :is_affine) && get_attribute(v, :is_affine)
    push!(properties_string, "affine")
  elseif has_attribute(v, :is_fano) && get_attribute(v, :is_fano)
    push!(properties_string, "Fano")
  elseif has_attribute(v, :is_projective) && get_attribute(v, :is_projective)
    push!(properties_string, "projective")
  elseif has_attribute(v, :is_complete) && get_attribute(v, :is_complete)
    push!(properties_string, "complete")
  end

  if has_attribute(v, :dim)
    push!(properties_string, string(dim(v)) * "-dimensional")
  end

  push!(properties_string, "toric variety")

  if has_attribute(v, :has_torusfactor)
    if get_attribute(v, :has_torusfactor)
      push!(properties_string, "with a torusfactor")
    else
      push!(properties_string, "without torusfactors")
    end
  end

  join(io, properties_string, " ")
end

Base.show(io::IO, ::MIME"text/plain", v::NormalToricVarietyType) = Base.show(pretty(io), v)
