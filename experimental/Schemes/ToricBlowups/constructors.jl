@doc raw"""
    blow_up(v::NormalToricVarietyType, new_ray::AbstractVector{<:IntegerUnion}; coordinate_name::String = "e")

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
Normal toric variety

julia> blow_down_morphism = blow_up(P3, [0, 1, 1])
Toric blowdown morphism

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
function blow_up(v::NormalToricVarietyType, new_ray::AbstractVector{<:IntegerUnion}; coordinate_name::String = "e")
  return ToricBlowdownMorphism(v, normal_toric_variety(star_subdivision(v, new_ray)), coordinate_name)
end

@doc raw"""
    blow_up(v::NormalToricVarietyType, n::Int; coordinate_name::String = "e")

Blow up the toric variety by subdividing the n-th cone in the list
of *all* cones of the fan of `v`. This cone need not be maximal.
This function returns the corresponding blowdown morphism.

By default, we pick "e" as the name of the homogeneous coordinate for
the exceptional divisor. As third optional argument one can supply
a custom variable name.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> blow_down_morphism = blow_up(P3, 5)
Toric blowdown morphism

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
function blow_up(v::NormalToricVarietyType, n::Int; coordinate_name::String = "e")
  return ToricBlowdownMorphism(v, normal_toric_variety(star_subdivision(v, n)), coordinate_name)
end

@doc raw"""
    blow_up(v::NormalToricVarietyType, I::MPolyIdeal; coordinate_name::String = "e")

Blow up the toric variety by subdividing the cone in the list
of *all* cones of the fan of `v` which corresponds to the
provided ideal `I`. Note that this cone need not be maximal.

By default, we pick "e" as the name of the homogeneous coordinate for
the exceptional divisor. As third optional argument one can supply
a custom variable name.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal toric variety

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

julia> I2 = ideal([x2 * x3])
ideal(x2*x3)

julia> b2P3 = blow_up(P3, I2);

julia> codomain(b2P3) == P3
true
```
"""
function blow_up(v::NormalToricVarietyType, I::MPolyIdeal; coordinate_name::String = "e")
    @req base_ring(I) == cox_ring(v) "The ideal must be contained in the cox ring of the toric variety"
    indices = [findfirst(y -> y == x, gens(cox_ring(v))) for x in gens(I)]
    if length(indices) == ngens(I) && !(nothing in indices)
      # We perform this blowup with toric techniques.
      rs = matrix(ZZ, rays(v))
      new_ray = vec(sum([rs[i,:] for i in indices]))
      new_ray = new_ray ./ gcd(new_ray)
      return ToricBlowdownMorphism(v, normal_toric_variety(star_subdivision(v, new_ray)), coordinate_name)
    else
      # We rely on advanced techniques to conduct this blowup (if available).
      return _generic_blow_up(v, I)
    end
end

function _generic_blow_up(v::Any, I::Any)
  error("Not yet supported")
end

function _generic_blow_up(v::NormalToricVarietyType, I::MPolyIdeal)
  return blow_up(ideal_sheaf(v, I))
end
