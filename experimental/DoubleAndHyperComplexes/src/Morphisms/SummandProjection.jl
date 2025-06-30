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
  check::Bool

  function SummandProjectionFactory(
      inc::AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType};
      check::Bool=true
    ) where {DomainType, CodomainType, MorphismType}
    return new{MorphismType}(inc, check)
  end
end

@attributes mutable struct SummandProjection{DomainType, CodomainType, MorphismType} <: AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType, SummandProjection{DomainType, CodomainType, MorphismType}}
  internal_morphism::HyperComplexMorphism{DomainType, CodomainType, MorphismType}

  function SummandProjection(
      inc::AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType};
      check::Bool=true
    ) where {DomainType, CodomainType, MorphismType}
    map_factory = SummandProjectionFactory(inc; check)
    
    internal_morphism = HyperComplexMorphism(
                                             codomain(inc), domain(inc), 
                                             map_factory, cached=true, 
                                             offset=-offset(inc)
                                            )
    return new{CodomainType, DomainType, MorphismType}(internal_morphism)
  end
end

function (fac::SummandProjectionFactory)(self::AbsHyperComplexMorphism, i::Tuple)
  inc = fac.inc[i]
  dom = domain(inc)::FreeMod
  cod = codomain(inc)::FreeMod
  img_gens = elem_type(dom)[]
  img_gens_inc = images_of_generators(inc)
  if fac.check
    is_injective(inc) || error("given map is not injective")
    A = sparse_matrix(inc)
    all(is_one(only(v)[2]) || is_one(-only(v)[2]) for v in A) || error("given map is not of the required form")
  end
  for g in gens(cod)
    j = findfirst(==(g), img_gens_inc)
    if isnothing(j)
      k = findfirst(==(-g), img_gens_inc)
      if isnothing(k)
        push!(img_gens, zero(dom))
      else
        push!(img_gens, -dom[k::Int])
      end
    else
      push!(img_gens, dom[j::Int])
    end
  end
  res = hom(cod, dom, img_gens)
  if fac.check && compose(inc, res) != identity_map(dom)
    error("given map failed to be of the required format.")
  end
  return res
end

function can_compute(fac::SummandProjectionFactory, self::AbsHyperComplexMorphism, i::Tuple)
  return can_compute_index(fac.inc, i)
end


underlying_morphism(phi::SummandProjection) = phi.internal_morphism
