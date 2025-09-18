function _find_blowup_coordinate_name(X::NormalToricVarietyType, coordinate_name::Union{VarName, Nothing} = nothing)
  if coordinate_name !== nothing
    internal_coordinate_name = Symbol(coordinate_name)
    @req !(internal_coordinate_name in coordinate_names(X)) "Coordinate name already exists"
    return internal_coordinate_name
  else
    return _find_blowup_coordinate_name(coordinate_names(X))
  end
end

function _find_blowup_coordinate_name(vs::Vector{Symbol})
  i = 1
  coordinate_name = :e
  while coordinate_name in vs
    coordinate_name = Symbol(:e, i)
    i = i+1
  end
  return coordinate_name
end


@doc raw"""
    blow_up(X::NormalToricVarietyType, r::AbstractVector{<:IntegerUnion}; coordinate_name::Union{VarName, Nothing} = nothing)

Star subdivide the fan $\Sigma$ of the toric variety $X$ along the
primitive vector `r` in the support of $\Sigma$ (Section
11.1 Star Subdivisions in [CLS11](@cite)).
This function returns the corresponding morphism.

By default, we pick "e" as the name of the homogeneous coordinate for
the exceptional prime divisor. As optional argument one can supply a
custom variable name.

# Examples

```jldoctest
julia> X = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> phi = blow_up(X, [0, 2, 3])
Toric blowup morphism

julia> Y = domain(phi)
Normal toric variety

julia> cox_ring(Y)
Multivariate polynomial ring in 5 variables over QQ graded by
  x1 -> [1 0]
  x2 -> [1 2]
  x3 -> [1 3]
  x4 -> [1 0]
  e -> [0 -1]
```
"""
function blow_up(X::NormalToricVarietyType, r::AbstractVector{<:IntegerUnion}; coordinate_name::Union{VarName, Nothing} = nothing)
  return ToricBlowupMorphism(X, r, _find_blowup_coordinate_name(X, coordinate_name))
end

@doc raw"""
    blow_up_along_minimal_supercone_coordinates(X::NormalToricVarietyType, minimal_supercone_coords::AbstractVector{<:RationalUnion}; coordinate_name::Union{VarName, Nothing} = nothing)

This method first constructs the primitive vector $r$ by calling
`standard_coordinates`, then blows up $X$ along $r$ using `blow_up`.

# Examples

```jldoctest
julia> X = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> phi = blow_up_along_minimal_supercone_coordinates(X, [2, 3, 0])
Toric blowup morphism
```
"""
function blow_up_along_minimal_supercone_coordinates(X::NormalToricVarietyType, minimal_supercone_coords::AbstractVector{<:RationalUnion}; coordinate_name::Union{VarName, Nothing} = nothing)
  coords = Vector{QQFieldElem}(minimal_supercone_coords)
  r_QQ = Vector{QQFieldElem}(standard_coordinates(polyhedral_fan(X), coords))
  r = primitive_generator(r_QQ)
  @req r == r_QQ "The argument `minimal_supercone_coords` must correspond to a primitive vector"
  phi = blow_up(X, r; coordinate_name)
  set_attribute!(phi, :minimal_supercone_coordinates_of_exceptional_ray, coords)
  return phi
end

@doc raw"""
    blow_up(X::NormalToricVarietyType, n::Int; coordinate_name::Union{VarName, Nothing} = nothing)

Blow up the toric variety $X$ with polyhedral fan $\Sigma$ by star
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
the exceptional prime divisor. As optional argument one can supply
a custom variable name.

# Examples
```jldoctest
julia> X = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> phi = blow_up(X, 5)
Toric blowup morphism

julia> Y = domain(phi)
Normal toric variety

julia> cox_ring(Y)
Multivariate polynomial ring in 5 variables over QQ graded by
  x1 -> [1 0]
  x2 -> [0 1]
  x3 -> [0 1]
  x4 -> [1 0]
  e -> [1 -1]
```
"""
function blow_up(X::NormalToricVarietyType, n::Int; coordinate_name::Union{VarName, Nothing} = nothing)
  coords = zeros(QQ, n_rays(X))
  for i in 1:number_of_rays(X)
    cones(X)[n, i] && (coords[i] = QQ(1))
  end
  exceptional_ray_scaled = standard_coordinates(polyhedral_fan(X), coords)
  exceptional_ray, scaling_factor = primitive_generator_with_scaling_factor(exceptional_ray_scaled)
  coords = scaling_factor * coords
  return blow_up_along_minimal_supercone_coordinates(X, coords; coordinate_name)
end
