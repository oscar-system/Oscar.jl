#= 
# A template for implementation of your own new morphism of hypercomplexes. 
#
# This is similar to the template in the `Objects` folder. See there for 
# further information and see the rest of the source for concrete implementation 
# examples.

struct MyNewMorphismFactory{MorphismType} <: HyperComplexMorphismFactory{MorphismType}
  # Fields needed for production

  function MyNewMorphismFactory(...)

    # Assuming that MorphismType has been extracted from the input
    return new{MorphismType}(...)
  end
end

function (fac::MyNewMorphismFactory)(self::AbsHyperComplexMorphism, i::Tuple)
  # Implement the production of the outgoing map
end

function can_compute(fac::MyNewMorphismFactory, self::AbsHyperComplexMorphism, i::Tuple)
  # Decide whether the outgoing map from index i can be computed
end


@attributes mutable struct MyNewMorphism{DomainType, CodomainType, MorphismType} <: AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType, MyNewMorphism{DomainType, CodomainType, MorphismType}}
  internal_morphism::HyperComplexMorphism{DomainType, CodomainType, MorphismType}
  # Further specialized fields
  # ...

  function MyNewMorphism(...)
    map_factory = MyNewMorphismFactory(...)

    # Assuming that the domain `dom` and the codomain `cod` have 
    # been extracted from the input
    internal_morphism = HyperComplexMorphism(dom, cod, map_factory, cached=true, offset=[0 for i in 1:dim(dom)])
    # Assuming that the types have been extracted from the input
    return new{DomainType, CodomainType, MorphismType}(internal_morphism, ...)
  end
end

underlying_morphism(phi::MyNewMorphism) = phi.internal_morphism
=#
