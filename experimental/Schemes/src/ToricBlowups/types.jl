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

  function ToricBlowupMorphism(X::NormalToricVarietyType, primitive_vector::AbstractVector{<:IntegerUnion}, coordinate_name::Symbol)
    # Construct the new variety
    Y = normal_toric_variety(star_subdivision(X, primitive_vector))

    # Compute position of the exceptional ray
    rays_Y = matrix(ZZ, rays(Y))
    index_of_exceptional_ray = findfirst(
      i -> primitive_vector == rays_Y[i, :],
      1:n_rays(Y),
    )
    @req index_of_exceptional_ray !== nothing "Could not identify position of exceptional ray"

    # Set variable names of Y
    var_names_X = symbols(cox_ring(X))
    @req !(coordinate_name in var_names_X) "The name for the blowup coordinate is already taken"
    var_names_Y = Vector{Symbol}(undef, n_rays(Y))
    rays_X = matrix(ZZ, rays(X))
    indices_X = Dict{AbstractVector, Int64}([rays_X[i,:]=>i for i in 1:n_rays(X)])
    for i in 1:n_rays(Y)
      if haskey(indices_X, rays_Y[i,:])
        var_names_Y[i] = var_names_X[indices_X[rays_Y[i,:]]]
      else
        var_names_Y[i] = coordinate_name
      end
    end
    set_attribute!(Y, :coordinate_names, var_names_Y)
    if n_rays(Y) > n_rays(X)
      @req coordinate_name in coordinate_names(Y) "Desired blowup variable name was not assigned"
    end

    # Construct the toric morphism
    phi_toric = toric_morphism(
      Y,
      identity_matrix(ZZ, ambient_dim(polyhedral_fan(X))),
      X;
      check=false,
    )

    # Construct the object
    phi = new{typeof(domain(phi_toric)), typeof(codomain(phi_toric))}(
      phi_toric, index_of_exceptional_ray
    )

    # Avoid recomputation
    if has_attribute(X, :has_torusfactor)
      set_attribute!(Y, :has_torusfactor, has_torusfactor(X))
    end
    if has_attribute(X, :is_orbifold) && is_orbifold(X)
      set_attribute!(Y, :is_orbifold, is_orbifold(X))
    end
    if has_attribute(X, :is_smooth) && is_smooth(X)
      if all(i -> primitive_vector[i] in [0, 1], 1:ambient_dim(X))
        set_attribute!(Y, :is_smooth, is_smooth(X))
      end
    end

    return phi
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
