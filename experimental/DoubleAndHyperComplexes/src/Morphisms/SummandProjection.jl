#= 
# Say you have an injective morphism of complexes `ι : C* ↪ D*` of `FreeMod`ules where 
#
#   ιᵏ : Cᵏ ≅ Rᵐ ↪ Rⁿ ≅ Dᵏ
#
# takes the generators of Cᵏ to a subset of the generators of Dᵏ for each k. 
# Then Rⁿ ≅ Rᵐ ⊕ Rⁿ⁻ᵐ decomposes as a direct sum. This type of morphism 
# computes the corresponding projection `π : D* → C*`. 
=#

struct SummandProjectionFactory{MorphismType} <: HyperComplexMorphismFactory{MorphismType}
  inc::AbsHyperComplexMorphism

  function SummandProjectionFactory(
      inc::AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType}
    ) where {DomainType, CodomainType, MorphismType}
    return new{MorphismType}(inc)
  end
end

function (fac::SummandProjectionFactory)(self::AbsHyperComplexMorphism, i::Tuple)
  inc = fac.inc[i]
  dom = domain(inc)::FreeMod
  cod = codomain(inc)::FreeMod
  img_gens = elem_type(dom)[]
  img_gens_inc = images_of_generators(inc)
  for (i, g) in enumerate(gens(cod))
    j = findfirst(==(g), img_gens_inc)
    if isnothing(j)
      push!(img_gens, zero(dom))
    else
      push!(img_gens, dom[j::Int])
    end
  end
  return hom(cod, dom, img_gens)
end

function can_compute(fac::SummandProjectionFactory, self::AbsHyperComplexMorphism, i::Tuple)
  return can_compute_index(fac.inc, i)
end


@attributes mutable struct SummandProjection{DomainType, CodomainType, MorphismType} <: AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType, SummandProjection{DomainType, CodomainType, MorphismType}}
  internal_morphism::HyperComplexMorphism{DomainType, CodomainType, MorphismType}

  function SummandProjection(
      inc::AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType}
    ) where {DomainType, CodomainType, MorphismType}
    map_factory = SummandProjectionFactory(inc)
    
    internal_morphism = HyperComplexMorphism(
                                             codomain(inc), domain(inc), 
                                             map_factory, cached=true, 
                                             offset=-offset(inc)
                                            )
    return new{CodomainType, DomainType, MorphismType}(internal_morphism)
  end
end

underlying_morphism(phi::SummandProjection) = phi.internal_morphism
