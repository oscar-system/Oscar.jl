@attributes mutable struct ToricBlowdownMorphism{
                                                 DomainType <: NormalToricVariety, 
                                                 CodomainType <: NormalToricVariety
    } <: AbsSimpleBlowdownMorphism{DomainType,
                                   CodomainType,
                                   Nothing,
                                   ToricBlowdownMorphism
                                  }
  toric_morphism::ToricMorphism
  new_ray::AbstractVector{<:IntegerUnion}

  # additional fields for caching
  center::IdealSheaf
  exceptional_divisor::ToricDivisor
  index_of_new_ray::Integer

  function ToricBlowdownMorphism(bl::ToricMorphism, new_ray::AbstractVector{<:IntegerUnion})
    return new{typeof(domain(bl)), typeof(codomain(bl))}(bl, new_ray)
  end

  function ToricBlowdownMorphism(bl::ToricMorphism, new_ray::AbstractVector{<:IntegerUnion}, center::IdealSheaf)
    return new{typeof(domain(bl)), typeof(codomain(bl))}(bl, new_ray, center)
  end

  function ToricBlowdownMorphism(
      Y::NormalToricVariety, new_ray::AbstractVector{<:IntegerUnion}, 
      coordinate_name::String, set_attributes::Bool
    )
    bl = _blow_up(Y, star_subdivision(Y, new_ray); 
                  coordinate_name = coordinate_name, 
                  set_attributes = set_attributes
                 )
    return new{typeof(domain(bl)), typeof(codomain(bl))}(bl, new_ray)
  end
end

ToricBlowdownMorphism(tbdm::ToricMorphism, new_ray::RayVector{QQFieldElem}) = ToricBlowdownMorphism(tbdm, _to_integer_vector(new_ray))

@doc raw"""
    blow_up(v::NormalToricVarietyType, new_ray::AbstractVector{<:IntegerUnion}; coordinate_name::String = "e", set_attributes::Bool = true)

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
function blow_up(v::NormalToricVarietyType, new_ray::AbstractVector{<:IntegerUnion}; coordinate_name::String = "e", set_attributes::Bool = true)
  return ToricBlowdownMorphism(v, new_ray, coordinate_name, set_attributes)
end

@doc raw"""
    blow_up(v::NormalToricVarietyType, n::Int; coordinate_name::String = "e", set_attributes::Bool = true)

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
function blow_up(v::NormalToricVarietyType, n::Int; coordinate_name::String = "e", set_attributes::Bool = true)
  # We assemble the center of the blowup
  c = cones(v)[n, :]::Polymake.SparseVectorBool
  g = gens(cox_ring(v))[c]
  I = ideal(cox_ring(v), g)
  II = IdealSheaf(v, I)
  
  bl = _blow_up(v, star_subdivision(v, n); coordinate_name = coordinate_name, set_attributes = set_attributes)
  new_rays = [v for v in rays(domain(bl)) if !(v in rays(codomain(bl)))]
  @assert length(new_rays) == 1 "there must be exactly one new ray"
  new_ray = first(new_rays)
  return ToricBlowdownMorphism(bl, _to_integer_vector(new_ray), II)
end

@doc raw"""
    blow_up(v::NormalToricVarietyType, I::MPolyIdeal; coordinate_name::String = "e", set_attributes::Bool = true)

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

julia> I2 = ideal([x2 * x3])
ideal(x2*x3)

julia> b2P3 = blow_up(P3, I2);

julia> codomain(b2P3) == P3
true
```
"""
function blow_up(v::NormalToricVarietyType, I::MPolyIdeal; coordinate_name::String = "e", set_attributes::Bool = true)
    @req base_ring(I) == cox_ring(v) "The ideal must be contained in the cox ring of the toric variety"
    indices = [findfirst(y -> y == x, gens(cox_ring(v))) for x in gens(I)]
    if length(indices) == ngens(I) && !(nothing in indices)
      # We perform this blowup with toric techniques.
      rs = matrix(ZZ, rays(v))
      new_ray = vec(sum([rs[i,:] for i in indices]))
      new_ray = new_ray ./ gcd(new_ray)
      II = IdealSheaf(v, I)
      bl = blow_up(v, new_ray; coordinate_name = coordinate_name, set_attributes = set_attributes)
      bl.center = II
      return bl
    else
      # We rely on advanced techniques to conduct this blowup (if available).
      return _generic_blow_up(v, I)
    end
end

function _generic_blow_up(v::Any, I::Any)
  error("Not yet supported")
end

function _blow_up(v::NormalToricVarietyType, new_fan::PolyhedralFan{QQFieldElem}; coordinate_name::String = "e", set_attributes::Bool = true)
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



# Forwarding all the essential functionality for morphisms of schemes 
underlying_morphism(bl::ToricBlowdownMorphism) = bl.toric_morphism
# now `domain`, `codomain`, `covering_morphism`, etc. should run automatically 
# due to our previous work.

new_ray(bl::ToricBlowdownMorphism) = bl.new_ray

function index_of_new_ray(bl::ToricBlowdownMorphism)
  if !isdefined(bl, :index_of_new_ray)
    X = domain(bl)
    v = new_ray(bl)
    j = findfirst(w->w == v, rays(X))
    @assert j!==nothing "ray not found"
    bl.index_of_new_ray = j
  end
  return bl.index_of_new_ray
end


# Implement the interface for AbsSimpleBlowdownMorphism
function center(bl::ToricBlowdownMorphism)
  if !isdefined(bl, :center)
    # TODO: The implementation below is highly inefficient. Improve it if you know how.
    X = domain(bl)
    S = cox_ring(X)
    x = gens(S)
    r = new_ray(bl)
    j = index_of_new_ray(bl)
    I = ideal(S, x[j])
    II = IdealSheaf(X, I)
    JJ = pushforward(bl, II)::IdealSheaf
    bl.center = JJ
  end
  return bl.center
end

function exceptional_divisor(bl::ToricBlowdownMorphism)
  if !isdefined(bl, :exceptional_divisor)
    X = domain(bl)
    S = cox_ring(X)
    x = gens(S)
    j = index_of_new_ray(bl)
    help_list = [i == j ? 1 : 0 for i in 1:ngens(S)]
    td = toric_divisor(X, help_list)
    @assert is_cartier(td) "exceptional divisor must be Cartier"
    @assert is_prime(td) "exceptional divisor must be prime"
    bl.exceptional_divisor = td
  end
  return bl.exceptional_divisor
end

function _to_vector(v::GrpAbFinGenElem)
  return [v[i] for i in 1:ngens(parent(v))]
end

function _to_integer_vector(v::RayVector{QQFieldElem})
  @assert all(x->isone(denominator(x)), v) "all denominators must be trivial"
  return [numerator(v[i]) for i in 1:length(v)]
end

########################################################################
# Forwarding of toric attributes                                       #
########################################################################
grid_morphism(bl::ToricBlowdownMorphism) = grid_morphism(underlying_morphism(bl))
morphism_on_torusinvariant_weil_divisor_group(bl::ToricBlowdownMorphism) = morphism_on_torusinvariant_weil_divisor_group(underlying_morphism(bl))
morphism_on_torusinvariant_cartier_divisor_group(bl::ToricBlowdownMorphism) = morphism_on_torusinvariant_cartier_divisor_group(underlying_morphism(bl))
morphism_on_class_group(bl::ToricBlowdownMorphism) = morphism_on_class_group(underlying_morphism(bl))
morphism_on_picard_group(bl::ToricBlowdownMorphism) = morphism_on_picard_group(underlying_morphism(bl))

#= For the future....
########################################################################
# Enabling the HasToricSubObjectTrait for ToricBlowdownMorphism        #
########################################################################
HasToricSubObjectTrait(::Type{T}) where {T<:ToricBlowdownMorphism} = HasToricSubObjectTrait{ToricMorphism}()

toric_sub_object(bl::ToricBlowdownMorphism) = underlying_morphism(bl)

########################################################################
# Enabling the functionality in general                                #
########################################################################
function grid_morphism(phi::Any)
  return _grid_morphism(HasToricSubObjectTrait(phi), phi)
end

function _grid_morphism(::HasToricSubObject, phi)
  return grid_morphism(toric_sub_object(phi))
end
=#


########################################################################
# Forwarding other functionality of the toric morphism                 #
########################################################################

# TODO: Arithmetic of toric morphisms.

########################################################################
# Further functionality required for the generic code of blowups
########################################################################
@attr IdealSheaf function ideal_sheaf(td::ToricDivisor)
  @assert is_cartier(td) "ideal sheaf can only be generated if the divisor is cartier"
  X = toric_variety(td)
  if is_prime(td)
    S = cox_ring(X)
    x = gens(S)
    j = findfirst(x->x==1, coefficients(td)) # Find out which one of the rays we actually have
    @assert j !== nothing "no ray was found"
    II = IdealSheaf(X, ideal(S, x[j]))
    return II
  end
  coeffs = coefficients(td)::Vector{ZZRingElem}
  prime_divisors = torusinvariant_prime_divisors(X)
  return prod(II^k for (II, k) in zip(ideal_sheaf.(prime_divisors), coeffs))
end



########################################################################
# Total transform
########################################################################

function total_transform(f::AbsSimpleBlowdownMorphism, II::IdealSheaf)
  return pullback(f, II)
end

function total_transform(f::AbsBlowdownMorphism, II::IdealSheaf)
  return pullback(f, II)
end
