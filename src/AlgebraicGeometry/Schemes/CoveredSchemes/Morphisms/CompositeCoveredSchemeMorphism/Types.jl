########################################################################
# Composite morphism of covered schemes
########################################################################

@doc raw"""
    CompositeCoveredSchemeMorphism{
        DomainType<:AbsCoveredScheme,
        CodomainType<:AbsCoveredScheme,
        BaseMorphismType
       } <: AbsCoveredSchemeMorphism{
                                 DomainType,
                                 CodomainType,
                                 BaseMorphismType,
                                 CoveredSchemeMorphism
                                }

A special concrete type of an `AbsCoveredSchemeMorphism` of the
form ``f = hᵣ ∘ hᵣ₋₁ ∘ … ∘ h₁: X → Y`` for arbitrary
`AbsCoveredSchemeMorphism`s ``h₁ : X → Z₁``, ``h₂ : Z₁ → Z₂``, ...,
``hᵣ : Zᵣ₋₁ → Y``.

Since every such morphism ``hⱼ`` will in general have an underlying
`CoveringMorphism` with `domain` and `codomain` `covering` actual
composition of such a sequence of morphisms will lead to an exponential
increase in complexity of these coverings because of the necessary
refinements. Nevertheless, the pullback or pushforward of various objects
on either ``X`` or ``Y`` through such a chain of maps is possible stepwise.
This type allows one to have one concrete morphism rather than a list
of morphisms and to reroute such calculations to iteration over the
various maps.

In addition to the usual functionality of the `AbsCoveredSchemeMorphism`
interface, this concrete type has the getters

    maps(f::CompositeCoveredSchemeMorphism)

to obtain a list of the ``hⱼ`` and `map(f, j)` to obtain the `j`-th map
directly.
"""
@attributes mutable struct CompositeCoveredSchemeMorphism{
    DomainType<:AbsCoveredScheme,
    CodomainType<:AbsCoveredScheme,
    BaseMorphismType
   } <: AbsCoveredSchemeMorphism{
                                 DomainType,
                                 CodomainType,
                                 BaseMorphismType,
                                 CoveredSchemeMorphism
                                }
  maps::Vector{<:AbsCoveredSchemeMorphism}

  # fields for caching
  composed_map::AbsCoveredSchemeMorphism

  function CompositeCoveredSchemeMorphism(maps::Vector{<:AbsCoveredSchemeMorphism})
    n = length(maps)
    for i in 1:n-1
      @assert codomain(maps[i]) === domain(maps[i+1]) "maps are not compatible"
    end
    # TODO: Take care of non-trivial base changes!
    return new{typeof(domain(first(maps))), typeof(codomain(maps[end])), Nothing}(maps)
  end
end
