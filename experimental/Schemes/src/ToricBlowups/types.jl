########################################################################
# Type definition                                                      #
########################################################################

@attributes mutable struct ToricBlowupMorphism{
  DomainType <: NormalToricVarietyType, 
  CodomainType <: NormalToricVarietyType,
  CenterDataType <: Union{
    AbstractVector{<:IntegerUnion},
    MPolyIdeal,
    ToricIdealSheafFromCoxRingIdeal,
    IdealSheaf,
  },
  CenterUnnormalizedType <: Union{
    ToricIdealSheafFromCoxRingIdeal,
    IdealSheaf,
  },
} <: AbsSimpleBlowupMorphism{DomainType, CodomainType, ToricBlowupMorphism}

  toric_morphism::ToricMorphism
  index_of_new_ray::Integer
  center_data::CenterDataType
  center_unnormalized::CenterUnnormalizedType
  exceptional_prime_divisor::ToricDivisor

  function _toric_blowup_morphism(v::NormalToricVarietyType, new_variety::NormalToricVarietyType, coordinate_name::String, new_ray::AbstractVector{<:IntegerUnion}, center_data::CenterDataType) where CenterDataType <: Union{
    AbstractVector{<:IntegerUnion},
    MPolyIdeal,
    ToricIdealSheafFromCoxRingIdeal,
    IdealSheaf,
  }
    # Compute position of new ray
    new_rays = matrix(ZZ, rays(new_variety))
    position_new_ray = findfirst(i->new_ray==new_rays[i,:], 1:n_rays(new_variety))
    @req position_new_ray !== nothing "Could not identify position of new ray"

    # Set variable names of the new variety
    old_vars = string.(symbols(cox_ring(v)))
    @req !(coordinate_name in old_vars) "The name for the blowup coordinate is already taken"
    new_vars = Vector{String}(undef, n_rays(new_variety))
    old_rays = matrix(ZZ, rays(v))
    old_indices = Dict{AbstractVector, Int64}([old_rays[i,:]=>i for i in 1:n_rays(v)])
    for i in 1:n_rays(new_variety)
      if haskey(old_indices, new_rays[i,:])
        new_vars[i] = old_vars[old_indices[new_rays[i,:]]]
      else
        new_vars[i] = coordinate_name
      end
    end
    set_attribute!(new_variety, :coordinate_names, new_vars)
    if n_rays(new_variety) > n_rays(v)
      @assert coordinate_name in coordinate_names(new_variety) "Desired blowup variable name was not assigned"
    end

    # Construct the toric morphism and construct the object
    bl_toric = toric_morphism(new_variety, identity_matrix(ZZ, ambient_dim(polyhedral_fan(v))), v; check=false)
    return bl_toric, position_new_ray, center_data
  end
  
  function ToricBlowupMorphism(v::NormalToricVarietyType, new_variety::NormalToricVarietyType, coordinate_name::String, new_ray::AbstractVector{<:IntegerUnion}, center_data::CenterDataType, center_unnormalized::ToricIdealSheafFromCoxRingIdeal) where CenterDataType <: Union{
    AbstractVector{<:IntegerUnion},
    MPolyIdeal,
    ToricIdealSheafFromCoxRingIdeal,
    IdealSheaf,
  }
    bl_toric, position_new_ray, center_data = _toric_blowup_morphism(v, new_variety, coordinate_name, new_ray, center_data)
    bl = new{
      typeof(domain(bl_toric)),
      typeof(codomain(bl_toric)),
      typeof(center_data),
      typeof(center_unnormalized),
    }(bl_toric, position_new_ray, center_data, center_unnormalized)
    if has_attribute(v, :has_torusfactor)
      set_attribute!(bl, :has_torusfactor, has_torusfactor(v))
    end
    return bl
  end
  
  function ToricBlowupMorphism(v::NormalToricVarietyType, new_variety::NormalToricVarietyType, coordinate_name::String, new_ray::AbstractVector{<:IntegerUnion}, center_data::CenterDataType) where CenterDataType <: Union{
    AbstractVector{<:IntegerUnion},
    MPolyIdeal,
    ToricIdealSheafFromCoxRingIdeal,
    IdealSheaf,
  }
    bl_toric, position_new_ray, center_data = _toric_blowup_morphism(v, new_variety, coordinate_name, new_ray, center_data)
    bl = new{
      typeof(domain(bl_toric)),
      typeof(codomain(bl_toric)),
      typeof(center_data),
      IdealSheaf{typeof(v), AbsAffineScheme, Ideal, Map},
    }(bl_toric, position_new_ray, center_data)
    if has_attribute(v, :has_torusfactor)
      set_attribute!(bl, :has_torusfactor, has_torusfactor(v))
    end
    return bl
  end
end



########################################################################
# Arithmetic for toric blowup moprhism and toric morphisms           #
########################################################################

function Base.:+(tm1::ToricBlowupMorphism, tm2::ToricBlowupMorphism)
  @req domain(tm1) === domain(tm2) "The morphisms must have identical domains"
  @req codomain(tm1) === codomain(tm2) "The morphisms must have identical codomains"
  return toric_morphism(domain(tm1), grid_morphism(tm1) + grid_morphism(tm2), codomain(tm1))
end

function Base.:-(tm1::ToricBlowupMorphism, tm2::ToricBlowupMorphism)
  @req domain(tm1) === domain(tm2) "The morphisms must have identical domains"
  @req codomain(tm1) === codomain(tm2) "The morphisms must have identical codomains"
  return toric_morphism(domain(tm1), grid_morphism(tm1) - grid_morphism(tm2), codomain(tm1))
end

function Base.:*(c::T, tm::ToricBlowupMorphism) where T <: IntegerUnion
new_grid_morphism = hom(domain(grid_morphism(tm)), codomain(grid_morphism(tm)), c * matrix(grid_morphism(tm)))
return toric_morphism(domain(tm), new_grid_morphism, codomain(tm))
end

Base.:+(tm1::ToricBlowupMorphism, tm2::ToricMorphism) = underlying_morphism(tm1) + tm2
Base.:-(tm1::ToricBlowupMorphism, tm2::ToricMorphism) = underlying_morphism(tm1) - tm2
Base.:+(tm1::ToricMorphism, tm2::ToricBlowupMorphism) = tm1 + underlying_morphism(tm2)
Base.:-(tm1::ToricMorphism, tm2::ToricBlowupMorphism) = tm1 - underlying_morphism(tm2)



######################################################
# Composition of toric blowups and toric morphisms #
######################################################

function Base.:*(tm1::ToricBlowupMorphism, tm2::ToricBlowupMorphism)
  @req codomain(tm1) === domain(tm2) "The codomain of the first morphism must be identically the same as the domain of the second morphism"
  return toric_morphism(domain(tm1), grid_morphism(tm1) * grid_morphism(tm2), codomain(tm2))
end

Base.:*(tm1::ToricMorphism, tm2::ToricBlowupMorphism) = tm1 * underlying_morphism(tm2)
Base.:*(tm1::ToricBlowupMorphism, tm2::ToricMorphism) = underlying_morphism(tm1) * tm2



####################################################
# Equality and hash of toric blowups             #
####################################################

function Base.:(==)(tm1::ToricBlowupMorphism, tm2::ToricBlowupMorphism)
  return domain(tm1) == domain(tm2) && codomain(tm1) == codomain(tm2) && grid_morphism(tm1) == grid_morphism(tm2)
end

function Base.hash(tm::ToricBlowupMorphism, h::UInt)
  b = 0x1a66f927cae2d409 % UInt
  h = hash(domain(tm), h)
  h = hash(codomain(tm), h)
  h = hash(grid_morphism(tm), h)
  return xor(h, b)
end



######################
# Display            #
######################

Base.show(io::IO, tbdm::ToricBlowupMorphism) = print(io, "Toric blowup morphism")
Base.show(io::IO, ::MIME"text/plain", tbdm::ToricBlowupMorphism) = Base.show(pretty(io), tbdm)
