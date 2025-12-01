#=
# Suppose f : C* → D* is a morphism of 1-dimensional complexes and one would like 
# to have a pseudo inverse for the induced map in homology f : H^i(C*) → H^i(D*).
# Then instead of spelling out the homologies completely, it might be beneficial to 
# work simply on representatives in the chain complex. This is what the maps here do.
=#
struct HomologyPseudoInverseFactory{MorphismType} <: HyperComplexMorphismFactory{MorphismType}
  f::AbsHyperComplexMorphism
  function HomologyPseudoInverseFactory(f::AbsHyperComplexMorphism)
    return new{Map}(f)
  end
end

function (fac::HomologyPseudoInverseFactory)(self::AbsHyperComplexMorphism, i::Tuple)
  J = collect(i) - offset(self)
  j = Tuple(J...)
  phi = fac.f[j]
  dom = domain(self) # codomain of fac.f
  F = dom[i]
  @assert codomain(phi) === F
  B, inc_B = boundary(dom, 1, i)
  Q, pr = quo(F, B)
  phi_red = compose(phi, pr)
  return map_from_func(F, domain(phi), x->preimage(phi_red, pr(x)))
end

function can_compute(fac::HomologyPseudoInverseFactory, self::AbsHyperComplexMorphism, i::Tuple)
  J = collect(i) - offset(self)
  j = Tuple(J...)
  return can_compute_index(fac.f, j)
end


@attributes mutable struct HomologyPseudoInverse{DomainType, CodomainType, MorphismType} <: AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType, HomologyPseudoInverse{DomainType, CodomainType, MorphismType}}
  internal_morphism::HyperComplexMorphism{DomainType, CodomainType, MorphismType}

  function HomologyPseudoInverse(
      f::AbsHyperComplexMorphism{DT, CT}
    ) where {DT, CT}
    @assert dim(domain(f)) == dim(codomain(f)) == 1
    map_factory = HomologyPseudoInverseFactory(f)
    internal_morphism = HyperComplexMorphism(codomain(f), domain(f), map_factory, cached=true, offset=-offset(f))
    # Assuming that the types have been extracted from the input
    return new{CT, DT, Map}(internal_morphism)
  end
end

underlying_morphism(phi::HomologyPseudoInverse) = phi.internal_morphism

