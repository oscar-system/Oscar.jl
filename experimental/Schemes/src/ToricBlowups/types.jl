########################################################################
# Type definition                                                      #
########################################################################

@attributes mutable struct ToricBlowupMorphism{
  DomainType <: NormalToricVarietyType,
  CodomainType <: NormalToricVarietyType,
} <: AbsSimpleBlowupMorphism{DomainType, CodomainType, ToricBlowupMorphism}
  toric_morphism::ToricMorphism
  index_of_exceptional_ray::Integer
  exceptional_prime_divisor::ToricDivisor

  function ToricBlowupMorphism(v::NormalToricVarietyType, exceptional_ray::AbstractVector{<:IntegerUnion}, coordinate_name::String)
    # Construct the new variety
    new_variety = normal_toric_variety(star_subdivision(v, exceptional_ray))

    # Compute position of the exceptional ray
    new_rays = matrix(ZZ, rays(new_variety))
    position_exceptional_ray = findfirst(
      i->exceptional_ray==new_rays[i,:],
      1:n_rays(new_variety)
    )
    @req position_exceptional_ray !== nothing "Could not identify position of exceptional ray"

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

    # Construct the toric morphism
    bl_toric = toric_morphism(
      new_variety,
      identity_matrix(ZZ, ambient_dim(polyhedral_fan(v))),
      v;
      check=false,
    )

    # Construct the object
    bl = new{typeof(domain(bl_toric)), typeof(codomain(bl_toric))}(
      bl_toric, position_exceptional_ray
    )

    # Avoid recomputation
    if has_attribute(v, :has_torusfactor)
      set_attribute!(bl, :has_torusfactor, has_torusfactor(v))
    end
    if has_attribute(v, :is_orbifold)
      set_attribute!(bl, :is_orbifold, is_orbifold(v))
    end
    if has_attribute(v, :is_smooth)
      if all(
        i -> rays(new_variety)[position_exceptional_ray][i] in [0, 1],
        1:ambient_dim(v),
      )
        set_attribute!(bl, :is_smooth, is_orbifold(v))
      end
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
  return toric_morphism(domain(tm1), lattice_homomorphism(tm1) + lattice_homomorphism(tm2), codomain(tm1))
end

function Base.:-(tm1::ToricBlowupMorphism, tm2::ToricBlowupMorphism)
  @req domain(tm1) === domain(tm2) "The morphisms must have identical domains"
  @req codomain(tm1) === codomain(tm2) "The morphisms must have identical codomains"
  return toric_morphism(domain(tm1), lattice_homomorphism(tm1) - lattice_homomorphism(tm2), codomain(tm1))
end

function Base.:*(c::T, tm::ToricBlowupMorphism) where T <: IntegerUnion
new_lattice_homomorphism = hom(domain(lattice_homomorphism(tm)), codomain(lattice_homomorphism(tm)), c * matrix(lattice_homomorphism(tm)))
return toric_morphism(domain(tm), new_lattice_homomorphism, codomain(tm))
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
  return toric_morphism(domain(tm1), lattice_homomorphism(tm1) * lattice_homomorphism(tm2), codomain(tm2))
end

Base.:*(tm1::ToricMorphism, tm2::ToricBlowupMorphism) = tm1 * underlying_morphism(tm2)
Base.:*(tm1::ToricBlowupMorphism, tm2::ToricMorphism) = underlying_morphism(tm1) * tm2



####################################################
# Equality and hash of toric blowups             #
####################################################

function Base.:(==)(tm1::ToricBlowupMorphism, tm2::ToricBlowupMorphism)
  return domain(tm1) == domain(tm2) && codomain(tm1) == codomain(tm2) && lattice_homomorphism(tm1) == lattice_homomorphism(tm2)
end

function Base.hash(tm::ToricBlowupMorphism, h::UInt)
  b = 0x1a66f927cae2d409 % UInt
  h = hash(domain(tm), h)
  h = hash(codomain(tm), h)
  h = hash(lattice_homomorphism(tm), h)
  return xor(h, b)
end



######################
# Display            #
######################

Base.show(io::IO, tbdm::ToricBlowupMorphism) = print(io, "Toric blowup morphism")
Base.show(io::IO, ::MIME"text/plain", tbdm::ToricBlowupMorphism) = Base.show(pretty(io), tbdm)
