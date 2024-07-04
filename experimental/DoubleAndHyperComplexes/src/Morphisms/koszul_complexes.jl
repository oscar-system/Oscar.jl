###
# Let ϕ : F → G be a morphism of free R-modules F ≅ Rᵐ and G ≅ Rⁿ. 
# Suppose that we have two elements v ∈ F and w ∈ G which give rise 
# to Koszul complexes 
#
#  0 → ⋀⁰ F → ⋀¹ F → ⋀² F → … → ⋀ᵐ⁻¹ F → ⋀ᵐ F → 0
#
# and
#
#  0 → ⋀⁰ G → ⋀¹ G → ⋀² G → … → ⋀ⁿ⁻¹ G → ⋀ⁿ G → 0
#
# given by exterior multiplication with v and w, respectively. 
# Then, in case ϕ(v) = w, there is an induced morphism of complexes 
# ⋀* F → ⋀* G coming from the maps on the exterior powers. 
#
# The struct below implements this morphism of complexes. 

struct KoszulKomplexMorphismFactory{MorphismType} <: HyperComplexMorphismFactory{MorphismType}
  original_map::MorphismType

  function KoszulKomplexMorphismFactory(
      phi::FreeModuleHom{DT, CT, Nothing}
    ) where {DT <: FreeMod, CT<:FreeMod}
    return new{typeof(phi)}(phi)
  end
end

function (fac::KoszulKomplexMorphismFactory)(self::AbsHyperComplexMorphism, i::Tuple)
  d = ngens(domain(fac.original_map)) - first(i)
  d == 0 && return hom(domain(self)[i], 
                       codomain(self)[ngens(codomain(fac.original_map)) - d],
                       [first(gens(codomain(self)[ngens(codomain(fac.original_map)) - d]))],
                       check=false
                      )
  return induced_map_on_exterior_power(fac.original_map, d, 
                                       domain=domain(self)[i], 
                                       codomain=codomain(self)[ngens(codomain(fac.original_map)) - d]
                                      )
end

function can_compute(fac::KoszulKomplexMorphismFactory, self::AbsHyperComplexMorphism, i::Tuple)
  first(i) <= ngens(domain(fac.original_map)) || return false
  ngens(domain(fac.original_map)) - first(i) <= ngens(codomain(fac.original_map)) || return false
  return true
end


@attributes mutable struct KoszulKomplexMorphism{DomainType, CodomainType, MorphismType} <: AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType, KoszulKomplexMorphism{DomainType, CodomainType, MorphismType}}
  internal_morphism::HyperComplexMorphism{DomainType, CodomainType, MorphismType}
  original_map::MorphismType

  function KoszulKomplexMorphism(
      dom::KoszulComplex, cod::KoszulComplex, phi::MorphismType
    ) where {MorphismType <: FreeModuleHom{<:FreeMod, <:FreeMod, Nothing}}
    v = defining_element(dom)
    w = defining_element(cod)
    F = parent(v)
    G = parent(w)
    @assert F === domain(phi)
    @assert G === codomain(phi)
    @assert phi(v) == w
    map_factory = KoszulKomplexMorphismFactory(phi)

    d = ngens(F) - ngens(G)
    internal_morphism = HyperComplexMorphism(dom, cod, 
                                             map_factory, cached=true, 
                                             offset=[d]
                                            )
    return new{typeof(dom), typeof(cod), MorphismType}(internal_morphism, phi)
  end
end

underlying_morphism(phi::KoszulKomplexMorphism) = phi.internal_morphism

function koszul_complex(
    phi::FreeModuleHom{DT, CT, Nothing},
    dom::KoszulComplex,
    cod::KoszulComplex
  ) where {DT<:FreeMod, CT<:FreeMod}
  return KoszulKomplexMorphism(dom, cod, phi)
end

  

