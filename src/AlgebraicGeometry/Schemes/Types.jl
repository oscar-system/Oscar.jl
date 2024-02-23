

### Abstract type for arbitrary schemes ###############################
#@doc raw"""
#    Scheme{BaseRingType<:Ring}
#
#A scheme over a ring ``𝕜`` of type `BaseRingType`.
#"""
#abstract type Scheme{BaseRingType} end
#
# Moved to src/forward_declarations.jl

@attr AffineScheme{S,S} function base_scheme(X::Scheme{S}) where {S<:Ring}
  return spec(base_ring(X))
end

### Abstract type for morphisms of arbitrary schemes ##################
@doc raw"""
    SchemeMor{DomainType, CodomainType, MorphismType, BaseMorType}

A morphism of schemes ``f : X → Y`` of type `MorphismType` with
``X`` of type `DomainType` and ``Y`` of type `CodomainType`.

When ``X`` and ``Y`` are defined over schemes ``BX`` and ``BY`` other
than ``Spec(𝕜)``, `BaseMorType` is the type of the underlying
morphism ``BX → BY``; otherwise, it can be set to `Nothing`.
"""
abstract type SchemeMor{
                        DomainType,
                        CodomainType,
                        MorphismType,
                        BaseMorType
                       } <: Map{
                                DomainType,
                                CodomainType,
                                SetMap,
                                MorphismType
                                }
end


### The empty scheme over a base ring #################################
struct EmptyScheme{BaseRingType}<:Scheme{BaseRingType}
  k::BaseRingType
  function EmptyScheme(k::BaseRingType) where {BaseRingType<:Ring}
    return new{BaseRingType}(k)
  end
end



########################################################################
# Abstract projective schemes                                          #
########################################################################
abstract type AbsProjectiveScheme{BaseRingType, RingType} <: Scheme{BaseRingType} end
