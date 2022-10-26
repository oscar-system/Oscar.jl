###########################################
# (1) General schemes
###########################################

export Scheme, SchemeMor, EmptyScheme

### Abstract type for arbitrary schemes ###############################
@Markdown.doc """
    Scheme{BaseRingType<:Ring} 

A scheme over a ring ``𝕜`` of type `BaseRingType`.
"""
abstract type Scheme{BaseRingType} end

### Abstract type for morphisms of arbitrary schemes ##################
@Markdown.doc """
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
                       } <: Hecke.Map{
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



###########################################
# (2) Types for affine schemes
###########################################

export AbsSpec, Spec, StdSpec

@Markdown.doc """
    AbsSpec{BaseRingType, RingType<:Ring}

An affine scheme ``X = Spec(R)`` with ``R`` of type `RingType` over
a ring ``𝕜`` of type `BaseRingType`.
"""
abstract type AbsSpec{BaseRingType, RingType<:Ring} <: Scheme{BaseRingType} end


@Markdown.doc """
    Spec{BaseRingType, RingType}

An affine scheme ``X = Spec(R)`` with ``R`` a Noetherian ring of type `RingType`
over a base ring ``𝕜`` of type `BaseRingType`.
"""
@attributes mutable struct Spec{BaseRingType, RingType} <: AbsSpec{BaseRingType, RingType}
  # the basic fields
  OO::RingType
  kk::BaseRingType

  function Spec(OO::MPolyQuoLocalizedRing)
    kk = coefficient_ring(base_ring(OO))
    return new{typeof(kk), typeof(OO)}(OO, kk)
  end
  function Spec(OO::MPolyLocalizedRing)
    kk = coefficient_ring(base_ring(OO))
    return new{typeof(kk), typeof(OO)}(OO, kk)
  end
  function Spec(OO::MPolyRing)
    kk = coefficient_ring(OO)
    return new{typeof(kk), typeof(OO)}(OO, kk)
  end
  function Spec(OO::MPolyQuo)
    kk = coefficient_ring(base_ring(OO))
    return new{typeof(kk), typeof(OO)}(OO, kk)
  end

  function Spec(R::Ring)
    return new{typeof(ZZ), typeof(R)}(R, ZZ)
  end

  function Spec(kk::Ring, R::Ring)
    return new{typeof(kk), typeof(R)}(R, kk)
  end

  function Spec(kk::Field)
    return new{typeof(kk), typeof(kk)}(kk, kk)
  end
end


StdSpec = AbsSpec{<:Ring, <:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}



########################################################################
# (3) Types for morphisms of affine schemes                            #
########################################################################

export AbsSpecMor, SpecMor
export OpenInclusion


@Markdown.doc """
    AbsSpecMor{DomainType<:AbsSpec,
               CodomainType<:AbsSpec,
               PullbackType<:Hecke.Map,
               MorphismType,
               BaseMorType
               }

Abstract type for morphisms ``f : X → Y`` of affine schemes where

  * ``X = Spec(S)`` is of type `DomainType`,
  * ``Y = Spec(R)`` is of type `CodomainType`,
  * ``f^* : R → S`` is a ring homomorphism of type `PullbackType`,
  * ``f`` itself is of type `MorphismType` (required for the Map interface),
  * if ``f`` is defined over a morphism of base schemes ``BX → BY``
    (e.g. a field extension), then this base scheme morphism is of
    type `BaseMorType`; otherwise, this can be set to `Nothing`.
"""
abstract type AbsSpecMor{
                         DomainType<:AbsSpec,
                         CodomainType<:AbsSpec,
                         PullbackType<:Hecke.Map,
                         MorphismType,
                         BaseMorType
                        }<:SchemeMor{DomainType, CodomainType, MorphismType, BaseMorType}
end


@Markdown.doc """
    SpecMor{DomainType<:AbsSpec,
            CodomainType<:AbsSpec,
            PullbackType<:Hecke.Map
           }

A morphism ``f : X → Y`` of affine schemes ``X = Spec(S)`` of type
`DomainType` and ``Y = Spec(R)`` of type `CodomainType`, both defined
over the same `base_ring`, with underlying ring homomorphism
``f^* : R → S`` of type `PullbackType`.
"""
@attributes mutable struct SpecMor{
                                   DomainType<:AbsSpec,
                                   CodomainType<:AbsSpec,
                                   PullbackType<:Hecke.Map
                                  } <: AbsSpecMor{DomainType,
                                                  CodomainType,
                                                  PullbackType,
                                                  SpecMor,
                                                  Nothing
                                                 }
  domain::DomainType
  codomain::CodomainType
  pullback::PullbackType

  function SpecMor(
      X::DomainType,
      Y::CodomainType,
      pullback::PullbackType;
      check::Bool=true
    ) where {DomainType<:AbsSpec, CodomainType<:AbsSpec, PullbackType<:Hecke.Map}
    OO(X) == codomain(pullback) || error("the coordinate ring of the domain does not coincide with the codomain of the pullback")
    OO(Y) == domain(pullback) || error("the coordinate ring of the codomain does not coincide with the domain of the pullback")
    if check
      # do some more expensive tests
    end
    return new{DomainType, CodomainType, PullbackType}(X, Y, pullback)
  end
end


@attributes mutable struct OpenInclusion{DomainType, CodomainType, PullbackType}<:AbsSpecMor{DomainType, CodomainType, PullbackType, OpenInclusion, Nothing}
  inc::SpecMor{DomainType, CodomainType, PullbackType}
  I::Ideal
  Z::Spec

  function OpenInclusion(f::AbsSpecMor, I::Ideal; check::Bool=true)
    U = domain(f)
    X = codomain(f)
    Z = subscheme(X, I)
    if check
      isempty(preimage(f, Z)) || error("image of the map is not contained in the complement of the vanishing locus of the ideal")
      #TODO: Do checks
    end
    return new{typeof(U), typeof(X), pullback_type(f)}(f, I, Z)
  end
end



########################################################################
# (4) Special type for OpenInclusion
########################################################################

export underlying_morphism

underlying_morphism(f::OpenInclusion) = f.inc
complement_ideal(f::OpenInclusion) = f.I
complement_scheme(f::OpenInclusion) = f.Z
