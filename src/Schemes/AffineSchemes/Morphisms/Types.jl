export AbsSpecMor, SpecMor
export OpenInclusion



########################################################################
# Abstract morphisms of affine schemes                                 #
########################################################################

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



########################################################################
# Minimal concrete type for morphisms of affine schemes                #
########################################################################

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

########################################################################
# A special type for open inclusions                                   #
########################################################################
@Markdown.doc """
    OpenInclusion{DomainType, CodomainType, PullbackType} <: AbsSpecMor

An open inclusion ``ι : U ↪ X`` of one affine scheme ``U`` into another 
one ``X``. 
"""
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


