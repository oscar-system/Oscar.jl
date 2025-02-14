function _find_blowup_coordinate_name(m::NormalToricVarietyType, coordinate_name::Union{String, Nothing} = nothing)
  if coordinate_name !== nothing
    @req !(coordinate_name in coordinate_names(m)) "Coordinate name already exists"
    return coordinate_name
  else
    return _find_blowup_coordinate_name(coordinate_names(m))
  end
end

function _find_blowup_coordinate_name(vs::Vector{String})
  i = 1
  coordinate_name = "e"
  while coordinate_name in vs
    coordinate_name = string("e", i)
    i = i+1
  end
  return coordinate_name
end


@doc raw"""
    blow_up(v::NormalToricVariety, exceptional_ray::AbstractVector{<:IntegerUnion}; coordinate_name::String)

Blow up the toric variety by subdividing the fan of the variety with the
provided exceptional ray. This function returns the corresponding morphism.

Note that this ray must be a primitive element in the lattice Z^d, with
d the dimension of the fan. In particular, it is currently impossible to
blow up along a ray which corresponds to a non-Q-Cartier divisor.

By default, we pick "e" as the name of the homogeneous coordinate for
the exceptional prime divisor. As optional argument one can supply a
custom variable name.

# Examples

```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> f = blow_up(P3, [0, 2, 3])
Toric blowup morphism

julia> bP3 = domain(f)
Normal toric variety

julia> cox_ring(bP3)
Multivariate polynomial ring in 5 variables over QQ graded by
  x1 -> [1 0]
  x2 -> [1 2]
  x3 -> [1 3]
  x4 -> [1 0]
  e -> [0 -1]
```
"""
function blow_up(v::NormalToricVarietyType, exceptional_ray::AbstractVector{<:IntegerUnion}; coordinate_name::Union{String, Nothing} = nothing)
  coordinate_name = _find_blowup_coordinate_name(v, coordinate_name)
  return ToricBlowupMorphism(v, exceptional_ray, coordinate_name)
end

@doc raw"""
    blow_up_along_minimal_supercone_coordinates(v::NormalToricVarietyType, minimal_supercone_coords::AbstractVector{<:IntegerUnion}; coordinate_name::Union{String, Nothing} = nothing)

This method first constructs the ray `r` by calling `standard_coordinates`, then blows up `v` along `r` using `blow_up`.

# Examples

```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> f = blow_up_along_minimal_supercone_coordinates(P2, [2, 3, 0])
Toric blowup morphism
```
"""
function blow_up_along_minimal_supercone_coordinates(v::NormalToricVarietyType, minimal_supercone_coords::AbstractVector{<:RationalUnion}; coordinate_name::Union{String, Nothing} = nothing)
  coords = Vector{QQFieldElem}(minimal_supercone_coords)
  ray_QQ = Vector{QQFieldElem}(standard_coordinates(polyhedral_fan(v), coords))
  ray = primitive_generator(ray_QQ)
  @assert ray == ray_QQ "The input vector must correspond to a primitive generator of a ray"
  phi = blow_up(v, ray; coordinate_name=coordinate_name)
  set_attribute!(phi, :minimal_supercone_coordinates_of_exceptional_ray, coords)
  return phi
end

@doc raw"""
    blow_up(v::NormalToricVariety, n::Int; coordinate_name::String = "e")

Blow up the toric variety $v$ with polyhedral fan $\Sigma$ by star
subdivision along the barycenter of the $n$-th cone $\sigma$ in the list
of all the cones of $\Sigma$.
We remind that the barycenter of a nonzero cone is the primitive
generator of the sum of the primitive generators of the extremal rays of
the cone (Exercise 11.1.10 in [CLS11](@cite)).
In the case all the cones of $\Sigma$ containing $\sigma$ are smooth,
this coincides with the star subdivision of $\Sigma$ relative to
$\sigma$ (Definition 3.3.17 of [CLS11](@cite)).
This function returns the corresponding morphism.

By default, we pick "e" as the name of the homogeneous coordinate for
the exceptional prime divisor. As third optional argument one can supply
a custom variable name.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> f = blow_up(P3, 5)
Toric blowup morphism

julia> bP3 = domain(f)
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
function blow_up(v::NormalToricVarietyType, n::Int; coordinate_name::Union{String, Nothing} = nothing)
  # minimal supercone coordinates
  coords = zeros(QQ, n_rays(v))
  for i in 1:number_of_rays(v)
    cones(v)[n, i] && (coords[i] = QQ(1))
  end
  exceptional_ray_scaled = standard_coordinates(polyhedral_fan(v), coords)
  exceptional_ray, scaling_factor = primitive_generator_with_scaling_factor(
    exceptional_ray_scaled
  )
  coords = scaling_factor * coords

  return blow_up_along_minimal_supercone_coordinates(v, coords; coordinate_name=coordinate_name)
end
