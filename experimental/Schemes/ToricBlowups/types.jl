########################################################################
# Type definition                                                      #
########################################################################

@attributes mutable struct ToricBlowdownMorphism{
  DomainType <: NormalToricVariety, 
  CodomainType <: NormalToricVariety} <: AbsSimpleBlowdownMorphism{DomainType, CodomainType, Nothing, ToricBlowdownMorphism}

  toric_morphism::ToricMorphism
  index_of_new_ray::Integer
  center::IdealSheaf
  exceptional_divisor::ToricDivisor

  function ToricBlowdownMorphism(v::NormalToricVariety, new_variety::NormalToricVariety, coordinate_name::String)
    old_vars = string.(symbols(cox_ring(v)))
    @req !(coordinate_name in old_vars) "The name for the blowup coordinate is already taken"
    new_vars = Vector{String}(undef, n_rays(v) + 1)
    index_of_new_ray = nothing
    for i in 1:n_rays(v)+1
        j = findfirst(==(rays(new_variety)[i]), rays(v))
        new_vars[i] = j !== nothing ? old_vars[j] : coordinate_name
        j === nothing && (index_of_new_ray = i)
    end
    @assert index_of_new_ray !==nothing "New ray not found -- some error in blow-up of toric variety"
    set_attribute!(new_variety, :coordinate_names, new_vars)
    bl = toric_morphism(new_variety, identity_matrix(ZZ, ambient_dim(polyhedral_fan(v))), v; check=false)
    return new{typeof(domain(bl)), typeof(codomain(bl))}(bl, index_of_new_ray)
  end
end

toric_blowdown_morphism(bl::ToricMorphism, new_ray::AbstractVector{<:IntegerUnion}, center::IdealSheaf) = ToricBlowdownMorphism(bl, new_ray, center)
toric_blowdown_morphism(Y::NormalToricVariety, new_ray::AbstractVector{<:IntegerUnion}, coordinate_name::String) = ToricBlowdownMorphism(Y, new_ray, coordinate_name)



########################################################################
# Arithmetic for toric blowdown moprhism and toric morphisms           #
########################################################################

function Base.:+(tm1::ToricBlowdownMorphism, tm2::ToricBlowdownMorphism)
  @req domain(tm1) === domain(tm2) "The morphisms must have identical domains"
  @req codomain(tm1) === codomain(tm2) "The morphisms must have identical codomains"
  return toric_morphism(domain(tm1), grid_morphism(tm1) + grid_morphism(tm2), codomain(tm1))
end

function Base.:-(tm1::ToricBlowdownMorphism, tm2::ToricBlowdownMorphism)
  @req domain(tm1) === domain(tm2) "The morphisms must have identical domains"
  @req codomain(tm1) === codomain(tm2) "The morphisms must have identical codomains"
  return toric_morphism(domain(tm1), grid_morphism(tm1) - grid_morphism(tm2), codomain(tm1))
end

function Base.:*(c::T, tm::ToricBlowdownMorphism) where T <: IntegerUnion
new_grid_morphism = hom(domain(grid_morphism(tm)), codomain(grid_morphism(tm)), c * matrix(grid_morphism(tm)))
return toric_morphism(domain(tm), new_grid_morphism, codomain(tm))
end

Base.:+(tm1::ToricBlowdownMorphism, tm2::ToricMorphism) = underlying_morphism(tm1) + tm2
Base.:-(tm1::ToricBlowdownMorphism, tm2::ToricMorphism) = underlying_morphism(tm1) - tm2
Base.:+(tm1::ToricMorphism, tm2::ToricBlowdownMorphism) = tm1 + underlying_morphism(tm2)
Base.:-(tm1::ToricMorphism, tm2::ToricBlowdownMorphism) = tm1 - underlying_morphism(tm2)



######################################################
# Composition of toric blowdowns and toric morphisms #
######################################################

function Base.:*(tm1::ToricBlowdownMorphism, tm2::ToricBlowdownMorphism)
  @req codomain(tm1) === domain(tm2) "The codomain of the first morphism must be identically the same as the domain of the second morphism"
  return toric_morphism(domain(tm1), grid_morphism(tm1) * grid_morphism(tm2), codomain(tm2))
end

Base.:*(tm1::ToricMorphism, tm2::ToricBlowdownMorphism) = tm1 * underlying_morphism(tm2)
Base.:*(tm1::ToricBlowdownMorphism, tm2::ToricMorphism) = underlying_morphism(tm1) * tm2



####################################################
# Equality and hash of toric blowdowns             #
####################################################

function Base.:(==)(tm1::ToricBlowdownMorphism, tm2::ToricBlowdownMorphism)
  return domain(tm1) == domain(tm2) && codomain(tm1) == codomain(tm2) && grid_morphism(tm1) == grid_morphism(tm2)
end

function Base.hash(tm::ToricBlowdownMorphism, h::UInt)
  b = 0x1a66f927cae2d409 % UInt
  h = hash(domain(tm), h)
  h = hash(codomain(tm), h)
  h = hash(grid_morphism(tm), h)
  return xor(h, b)
end



######################
# Display            #
######################

Base.show(io::IO, tbdm::ToricBlowdownMorphism) = print(io, "Toric blowdown morphism")
Base.show(io::IO, ::MIME"text/plain", tbdm::ToricBlowdownMorphism) = Base.show(pretty(io), tbdm)
