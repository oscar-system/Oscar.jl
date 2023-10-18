### Abstract type and interface for morphisms of complexes
abstract type AbsMorphismOfComplexes{DomainChainType, CodomainChainType, MorphismType} end

function domain(moc::AbsMorphismOfComplexes) 
  return domain(underlying_morphism(moc))
end

function codomain(moc::AbsMorphismOfComplexes) 
  return codomain(underlying_morphism(moc))
end

# Return the map 
function getindex(moc::AbsMorphismOfComplexes, i::Int)
  return underlying_morphism(moc)[i]
end

# Return the shift `k` such that the `i`-th morphism goes from C_i to D_{i + k}.
function degree(moc::AbsMorphismOfComplexes)
  return degree(underlying_morphism(moc))
end

# Return whether there is a known upper bound for the 
# range of maps for this morphism of complexes. 
#
# If the answer is `false`, but `extends_right` returns 
# true, then it must be safe to assume that every request
# for a map `moc[i]` is safe for i ≫ 0. 
#
# If, on the contrary, `has_upper_bound` is `true` and 
# `upper_bound` returns a number `k`, then one may 
# ask for `moc[i]` for `i ≤ k` and assume that the 
# maps beyond that bound are implicitly zero.
function has_upper_bound(moc::AbsMorphismOfComplexes)
  return has_upper_bound(underlying_morphism(moc))
end

function upper_bound(moc::AbsMorphismOfComplexes)
  return upper_bound(underlying_morphism(moc))
end

function has_lower_bound(moc::AbsMorphismOfComplexes)
  return has_lower_bound(underlying_morphism(moc))
end

function lower_bound(moc::AbsMorphismOfComplexes)
  return lower_bound(underlying_morphism(moc))
end

# Return whether this morphism of complexes knows how 
# to extend its range of maps to the left or the right.
function extends_left(moc::AbsMorphismOfComplexes)
  return extends_left(underlying_morphism(moc))
end

function extends_right(moc::AbsMorphismOfComplexes)
  return extends_right(underlying_morphism(moc))
end

function Base.range(moc::AbsMorphismOfComplexes)
  return range(underlying_morphism(moc))
end


### Abstract type and interface for the generation of morphisms

abstract type AbsMorphismFactory{MorphismType} end

# Construct the morphism from the i-th entry of the domain
function (fac::AbsMorphismFactory)(moc::AbsMorphismOfComplexes, i::Int)
  error("generation of morphisms not implemented for factories of type $(typeof(fac)); see the source code for details")
end


### Minimal concrete type for morphisms of complexes
mutable struct MorphismOfComplexes{DomainChainType, CodomainChainType, MorphismType} <: AbsMorphismOfComplexes{DomainChainType, CodomainChainType, MorphismType} 
  domain::ComplexOfMorphisms{DomainChainType}
  codomain::ComplexOfMorphisms{CodomainChainType}
  morphisms::Dict{Int, MorphismType}
  degree::Int

  # Production of the maps
  morphism_factory::AbsMorphismFactory{MorphismType}

  # Further attributes
  extends_right::Bool
  extends_left::Bool
  right_bound::Int
  left_bound::Int

  function MorphismOfComplexes(
      dom::ComplexOfMorphisms{T1}, cod::ComplexOfMorphisms{T2}, 
      morphism_factory::AbsMorphismFactory{MT}, degree::Int;
      right_bound::Union{Int, Nothing}=nothing, # By default we assume that the
      left_bound::Union{Int, Nothing}=nothing,  # morphisms are filled up with 
      extends_left::Bool=true,                  # zeroes on either side on request.
      extends_right::Bool=true                  # This depends on the factory used. 
    ) where {T1, T2, MT}
    mor_dict = Dict{Int, MT}()
    result = new{T1, T2, MT}(dom, cod, mor_dict, degree, morphism_factory)

    if right_bound !== nothing
      result.right_bound = right_bound
    end

    if left_bound !== nothing
      result.left_bound = left_bound
    end

    result.extends_left = extends_left
    result.extends_right= extends_right

    return result
  end
end

mutable struct LiftingMorphismsThroughResolution{T} <: AbsMorphismFactory{T}
   original_map::ModuleFPHom
  
   function LiftingMorphismsThroughResolution(phi::ModuleFPHom)
    return new{ModuleFPHom}(phi)
  end
end

