#= 
# Induced maps on tensor products.
#
# Say `ϕₖ: Cₖ* → Dₖ*` is a collection of `AbsHyperComplexMorphism`s for 
# `k = 1,…,n`, then we wish to compute the induced map
#
#   ϕ₁⊗ … ⊗ ϕₙ : C₁⊗ …⊗ Cₙ → D₁ ⊗ … ⊗ Dₙ
#
# on the tensor product complexes.
=#

struct InducedTensorProductMorphismFactory{MorphismType} <: HyperComplexMorphismFactory{MorphismType}
  factors::Vector{<:AbsHyperComplexMorphism}
  domain::HCTensorProductComplex
  codomain::HCTensorProductComplex

  function InducedTensorProductMorphismFactory(
      factors::Vector{<:AbsHyperComplexMorphism},
      domain::HCTensorProductComplex,
      codomain::HCTensorProductComplex
    )
    return new{ModuleFPHom}(factors, domain, codomain)
  end
end

function (fac::InducedTensorProductMorphismFactory)(self::AbsHyperComplexMorphism, i::Tuple)
  dims = [dim(c) for c in factors(fac.domain)]
  inds = Vector{Vector{Int}}()
  d0 = 0
  for d in dims
    push!(inds, [i[k] for k in d0+1:d0+d])
    d0 += d
  end
  maps = [fac.factors[i][inds[i]...] for i in 1:length(dims)]
  return tensor_product(maps; domain=fac.domain[i], codomain=fac.codomain[i])
end

function can_compute(fac::InducedTensorProductMorphismFactory, self::AbsHyperComplexMorphism, i::Tuple)
  return can_compute_index(fac.domain, i)
end


@attributes mutable struct InducedTensorProductMorphism{DomainType, CodomainType, MorphismType} <: AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType, InducedTensorProductMorphism{DomainType, CodomainType, MorphismType}}
  internal_morphism::HyperComplexMorphism{DomainType, CodomainType, MorphismType}

  function InducedTensorProductMorphism(
      factors::Vector{<:AbsHyperComplexMorphism};
      domain::HCTensorProductComplex = tensor_product([domain(phi) for phi in factors]),
      codomain::HCTensorProductComplex = tensor_product([codomain(phi) for phi in factors])
    )
    nfac = length(factors)
    @assert length(factors) == length(Oscar.factors(domain)) == length(Oscar.factors(codomain))
    @assert all(Oscar.domain(factors[i]) === Oscar.factors(domain)[i] for i in 1:nfac)
    @assert all(Oscar.codomain(factors[i]) === Oscar.factors(codomain)[i] for i in 1:nfac)
    map_factory = InducedTensorProductMorphismFactory(factors, domain, codomain)

    internal_morphism = HyperComplexMorphism(domain, codomain, map_factory, cached=true, offset=[0 for i in 1:dim(domain)])
    return new{typeof(domain), typeof(codomain), ModuleFPHom}(internal_morphism)
  end
end

underlying_morphism(phi::InducedTensorProductMorphism) = phi.internal_morphism

function tensor_product(
    factors::Vector{<:AbsHyperComplexMorphism};
    domain::HCTensorProductComplex = tensor_product([domain(phi) for phi in factors]),
    codomain::HCTensorProductComplex = tensor_product([codomain(phi) for phi in factors])
  )
  return InducedTensorProductMorphism(factors; domain, codomain)
end

