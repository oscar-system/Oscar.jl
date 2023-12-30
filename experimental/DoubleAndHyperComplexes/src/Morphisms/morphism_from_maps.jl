struct MorphismFromDictFactory{MorphismType} <: HyperComplexMorphismFactory{MorphismType}
  maps::Dict{<:Tuple, <:MorphismType}

  function MorphismFromDictFactory(
      ::Type{MorphismType}, 
      maps::Dict{<:Tuple, <:MorphismType}
    ) where {MorphismType}
    return new{MorphismType}(maps)
  end
end

function (fac::MorphismFromDictFactory)(self::AbsHyperComplexMorphism, i::Tuple)
  return fac.maps[i]
end

function can_compute(fac::MorphismFromDictFactory, self::AbsHyperComplexMorphism, i::Tuple)
  return haskey(fac.maps, i)
end


@attributes mutable struct MorphismFromDict{DomainType, CodomainType, MorphismType} <: AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType, MorphismFromDict{DomainType, CodomainType, MorphismType}}
  internal_morphism::HyperComplexMorphism{DomainType, CodomainType, MorphismType}

  function MorphismFromDict(
      dom::AbsHyperComplex{CTD, MTD},
      cod::AbsHyperComplex{CTC, MTC},
      maps::Dict{<:Tuple, MT};
      offset::Vector{Int} = [0 for i in 1:dim(dom)]
    ) where {CTD, MTD, CTC, MTC, MT}
    map_factory = MorphismFromDictFactory(MT, maps)
    @assert dim(dom) == dim(cod)
    @assert all(k->length(k) == dim(dom), keys(maps))
    for (k, phi) in maps
      @assert domain(phi) === dom[k]
      @assert codomain(phi) === cod[Tuple(collect(k) + offset)]
    end
    internal_morphism = HyperComplexMorphism(dom, cod, map_factory, cached=false, offset=offset)
    return new{typeof(dom), typeof(cod), MT}(internal_morphism)
  end
end

underlying_morphism(phi::MorphismFromDict) = phi.internal_morphism
