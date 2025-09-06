struct HomChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  d::AbsHyperComplex
  c::AbsHyperComplex
  auto_extend::Bool

  function HomChainFactory(
      ::Type{ChainType}, d::AbsHyperComplex, c::AbsHyperComplex;
      auto_extend::Bool=false
    ) where {ChainType}
    return new{ChainType}(d, c, auto_extend)
  end
end

function (fac::HomChainFactory)(hc::AbsHyperComplex, i::Tuple)
  I = collect(i)
  d = fac.d
  c = fac.c
  i1 = I[1:dim(d)]
  i2 = I[dim(d)+1:end]
  return hom(d[Tuple(-i1)], c[Tuple(i2)])[1]
end

function can_compute(fac::HomChainFactory, hc::AbsHyperComplex, i::Tuple)
  I = collect(i)
  d = fac.d
  c = fac.c
  i1 = I[1:dim(d)]
  i2 = I[dim(d)+1:end]
  return can_compute_index(fac.d, Tuple(-i1)) && can_compute_index(fac.c, Tuple(i2))
end


struct HomMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType}
  d::AbsHyperComplex
  c::AbsHyperComplex
  auto_extend::Bool

  function HomMapFactory(
      ::Type{MorphismType}, d::AbsHyperComplex, c::AbsHyperComplex;
      auto_extend::Bool=false
    ) where {MorphismType}
    return new{MorphismType}(d, c, auto_extend)
  end
end

function (fac::HomMapFactory)(hc::AbsHyperComplex, p::Int, i::Tuple)
  I = collect(i)
  d = fac.d
  c = fac.c
  i1 = I[1:dim(d)]
  i2 = I[dim(d)+1:end]

  dom = hc[i]
  inc = (direction(hc, p) == :chain ? -1 : 1)
  i_inc = Tuple(I + [k == p ? inc : 0 for k in 1:dim(hc)])
  cod = hc[i_inc]

  if iszero(dom) || iszero(cod)
    return hom(dom, cod, elem_type(cod)[zero(cod) for i in 1:ngens(dom)]; check=false)
  end

  if p <= dim(d)
    # contravariant induced map on first argument
    i1_inc = Tuple(-i1 - [k == p ? inc : 0 for k in 1:dim(d)])
    img_gens = [homomorphism_to_element(cod, compose(map(d, p, i1_inc), element_to_homomorphism(g))) for g in gens(dom)]
    return hom(dom, cod, img_gens; check=false)
  else
    # covariant induced map on second argument
    img_gens = [homomorphism_to_element(cod, compose(element_to_homomorphism(g), map(c, p - dim(d), Tuple(i2)))) for g in gens(dom)]
    return hom(dom, cod, img_gens; check=false)
  end
end

function can_compute(fac::HomMapFactory, h::AbsHyperComplex, p::Int, i::Tuple)
  fac.auto_extend && return true
  I = collect(i)
  d = fac.d
  c = fac.c

  # Check for computability of the domain and codomain
  can_compute_index(h, i) || return false
  J = I + (direction(h, p) == :chain ? -1 : 1)*[k==p ? 1 : 0 for k in 1:dim(h)]
  can_compute_index(h, Tuple(J)) || return false

  # Check for computability of the original map inducing the one here
  j1 = J[1:dim(d)]
  j2 = J[dim(d)+1:end]
  if p <= dim(fac.d)
    can_compute_map(d, p, Tuple(-j1)) || return false
  else
    can_compute_map(c, p-dim(d), Tuple(j2)) || return false
  end
  return true
end

### A concrete wrapper type for hom-complexes. 
# This allows for specialized functionality for this kind of complexes.
@attributes mutable struct HomComplex{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType} 
  domain::AbsHyperComplex
  codomain::AbsHyperComplex
  internal_complex::HyperComplex{ChainType, MorphismType}

  function HomComplex(
      dom::AbsHyperComplex, cod::AbsHyperComplex,
      internal_complex::HyperComplex{ChainType, MorphismType}) where {ChainType, MorphismType}
    return new{ChainType, MorphismType}(dom, cod, internal_complex)
  end
end

underlying_complex(c::HomComplex) = c.internal_complex

domain(c::HomComplex) = c.domain
codomain(c::HomComplex) = c.codomain

### User facing constructors
function hom(
    d::AbsHyperComplex{DCT, DMT}, c::AbsHyperComplex{CCT, CMT};
    auto_extend::Bool=true
  ) where {DCT, DMT, CCT, CMT}
  return _hom(d, c; auto_extend)
end

function _hom(
    d::AbsHyperComplex{DCT, DMT}, c::AbsHyperComplex{CCT, CMT};
    auto_extend::Bool=true
  ) where {DCT, DMT, CCT, CMT}
  NCT = typejoin(chain_type(d), chain_type(c))
  NMT = typejoin(morphism_type(d), morphism_type(c))
  
  d1 = dim(d)
  d2 = dim(c)

  chain_fac = HomChainFactory(NCT, d, c; auto_extend)
  map_fac = HomMapFactory(NMT, d, c; auto_extend)

  directions = [direction(d, i) for i in 1:d1]
  directions = vcat(directions, [direction(c, i) for i in 1:d2])

  upper_bounds = Vector{Union{Int, Nothing}}([(has_lower_bound(d, i) ? -lower_bound(d, i) : nothing) for i in 1:d1])
  upper_bounds = vcat(upper_bounds, Vector{Union{Int, Nothing}}([(has_upper_bound(c, i) ? upper_bound(c, i) : nothing) for i in 1:d2]))

  lower_bounds = Vector{Union{Int, Nothing}}([(has_upper_bound(d, i) ? -upper_bound(d, i) : nothing) for i in 1:d1])
  lower_bounds = vcat(lower_bounds, Vector{Union{Int, Nothing}}([(has_lower_bound(c, i) ? lower_bound(c, i) : nothing) for i in 1:d2]))


  return HomComplex(d, c,
                    HyperComplex(d1+d2, chain_fac, map_fac, directions, upper_bounds=upper_bounds, lower_bounds=lower_bounds)
                   )
end

#function hom(c::AbsSimpleComplex{ChainType}, M::ModuleFP; auto_extend::Bool=false) where {ChainType<:ModuleFP}
#  return SimpleComplexWrapper(_hom(c, ZeroDimensionalComplex(M); auto_extend))
#end

function hom(c::AbsHyperComplex{ChainType}, M::ModuleFP; auto_extend::Bool=false) where {ChainType<:ModuleFP}
  return _hom(c, ZeroDimensionalComplex(M); auto_extend)
end

#function hom(M::ModuleFP, c::AbsSimpleComplex{ChainType}; auto_extend::Bool=false) where {ChainType <: ModuleFP}
#  return SimpleComplexWrapper(_hom(ZeroDimensionalComplex(M), c; auto_extend))
#end

function hom(M::ModuleFP, c::AbsHyperComplex{ChainType}; auto_extend::Bool=false) where {ChainType <: ModuleFP}
  return _hom(ZeroDimensionalComplex(M), c; auto_extend)
end

#function hom(c1::AbsSimpleComplex, c2::AbsSimpleComplex; auto_extend::Bool=false)
#  return DoubleComplexWrapper(_hom(c1, c2; auto_extend))
#end

### Interpretation map for hom complexes
# An element `v` in the i-th homology of a total complex `tot` of a 
# hom complex `hom(a, b)` represents a homomorphism 
#
#   f : tot(a) -> tot(b)[i]
#
# from the total complex of `a` to the total complex of `b` shifted by i.
struct InterpretationMorphismFactory{MorphismType} <: HyperComplexMorphismFactory{MorphismType}
  tot::TotalComplex{<:ModuleFP, <:ModuleFPHom, <:HomComplex}
  dom::TotalComplex
  cod::TotalComplex
  v::ModuleFPElem
  i::Int

  function InterpretationMorphismFactory(
      tot::TotalComplex{ChainType, MorphismType, <:HomComplex}, 
      dom::TotalComplex, 
      cod::TotalComplex,
      v::ModuleFPElem, i::Int
    ) where {ChainType, MorphismType}
    @assert original_complex(dom) === domain(original_complex(tot))
    @assert original_complex(cod) === codomain(original_complex(tot))
    @assert parent(v) === tot[i]
    
    return new{MorphismType}(tot, dom, cod, v, i)
  end
end

function (fac::InterpretationMorphismFactory)(self::AbsHyperComplexMorphism, I::Tuple)
  i = first(I)
  tot = fac.tot
  dom = fac.dom
  cod = fac.cod
  offset = fac.i
  tot_orig = original_complex(tot)
  dom_orig = original_complex(dom)
  cod_orig = original_complex(cod)
  indices_dom = indices_in_summand(dom, i)
  j = i + offset
  indices_cod = indices_in_summand(cod, j)
  indices_tot = indices_in_summand(tot, offset)
  projections_dom = projections_for_summand(dom, i)
  injections_cod = injections_for_summand(cod, j)
  projections_tot = projections_for_summand(tot, offset)
  M = dom[i]
  N = cod[j]
  result = hom(M, N, elem_type(N)[zero(N) for i in 1:ngens(M)]; check=false)
  for (k, K) in enumerate(indices_dom)
    pr = projections_dom[k]
    for (l, L) in enumerate(indices_cod)
      KL = Tuple(vcat(-collect(K), collect(L)))
      kl = findfirst(==(KL), indices_tot)
      kl === nothing && error("index not found")
      v_loc = (projections_tot[kl])(fac.v)
      phi_loc = element_to_homomorphism(v_loc)
      @assert domain(phi_loc) === codomain(pr)
      inc = injections_cod[l]
      @assert codomain(phi_loc) === domain(inc)
      psi_loc = compose(compose(pr, phi_loc), inc)
      result = result + psi_loc
    end
  end

  # Our sign convention for total complexes makes the following adjustment 
  # necessary. See src/Objects/total_complexes.jl for details.
  return (-1)^(div(i + offset + dim(dom), 2)) * result
end

function can_compute(fac::InterpretationMorphismFactory, self::AbsHyperComplexMorphism, i::Tuple)
  can_compute_index(fac.dom, i) || return false
  can_compute_index(fac.cod, Tuple(i[1] + fac.i)) || return false
  can_compute_index(fac.tot, Tuple(fac.i)) || return false
  return true
end


@attributes mutable struct InterpretationMorphism{DomainType, CodomainType, MorphismType} <: AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType, InterpretationMorphism{DomainType, CodomainType, MorphismType}}
  internal_morphism::HyperComplexMorphism{DomainType, CodomainType, MorphismType}
  tot::TotalComplex{<:ModuleFP, <:ModuleFPHom, <:HomComplex}
  v::ModuleFPElem

  function InterpretationMorphism(
      tot::TotalComplex{ChainType, MorphismType, <:HomComplex}, 
      dom::TotalComplex, 
      cod::TotalComplex,
      v::ModuleFPElem, i::Int
    ) where {ChainType, MorphismType}
    map_factory = InterpretationMorphismFactory(tot, dom, cod, v, i)

    # Assuming that the domain `dom` and the codomain `cod` have 
    # been extracted from the input
    internal_morphism = HyperComplexMorphism(dom, cod, map_factory, cached=true, offset=[i])
    # Assuming that `MorphismType` has been extracted from the input
    return new{typeof(dom), typeof(cod), MorphismType}(internal_morphism, tot, v)
  end
end

underlying_morphism(phi::InterpretationMorphism) = phi.internal_morphism
module_element(phi::InterpretationMorphism) = phi.v

### The following function produces another function `interp` 
# which takes an element `v` of the `i`-th kernel module of `tot`
# and produces the associated homomorphism of (total-)complexes
#
#  f : tot(a) -> tot(b)
#
# where `tot = tot(hom(a, b))`.
function element_to_homomorphism_map(
    tot::TotalComplex{<:ModuleFP, <:ModuleFPHom, <:HomComplex}, 
    i::Int;
    domain::TotalComplex=total_complex(domain(original_complex(tot))),
    codomain::TotalComplex=total_complex(codomain(original_complex(tot)))
  )
  @assert original_complex(domain) === Oscar.domain(original_complex(tot))
  @assert original_complex(codomain) === Oscar.codomain(original_complex(tot))
  return function(v::SubquoModuleElem)
    parent(v) === kernel(tot, i)[1] || error("element does not belong to the correct module")
    w = ambient_representative(v)
    hom_comp = original_complex(tot)
    dom_tot = total_complex(domain)
    cod_tot = total_complex(codomain)
    return InterpretationMorphism(tot, domain, codomain, w, i)
  end
end

########################################################################
# Induced morphisms on hom-complexes                                   #
########################################################################
struct HomComplexMorphismFactory{MorphismType} <: HyperComplexMorphismFactory{MorphismType}
  domain_morphism::Union{Nothing, <:AbsHyperComplexMorphism}
  codomain_morphism::Union{Nothing, <:AbsHyperComplexMorphism}

  function HomComplexMorphismFactory(
      ::Type{MorphismType},
      dom_mor::Union{Nothing, <:AbsHyperComplexMorphism},
      cod_mor::Union{Nothing, <:AbsHyperComplexMorphism}
    ) where {MorphismType}
    return new{MorphismType}(dom_mor, cod_mor)
  end
end

function (fac::HomComplexMorphismFactory)(self::AbsHyperComplexMorphism, I::Tuple)
  if fac.domain_morphism === nothing && fac.codomain_morphism === nothing
    return id_hom(domain(self)[I])
  end

  i = collect(I)
  dom = domain(self)::HomComplex
  dom_dom = domain(dom)::AbsHyperComplex
  dom_cod = codomain(dom)::AbsHyperComplex
  cod = codomain(self)::HomComplex
  cod_dom = domain(cod)::AbsHyperComplex
  cod_cod = codomain(cod)::AbsHyperComplex
  d1 = fac.domain_morphism !== nothing ? offset(fac.domain_morphism) : [0 for i in 1:dim(dom_dom)]
  d2 = fac.codomain_morphism !== nothing ? offset(fac.codomain_morphism) : [0 for i in 1:dim(cod_cod)]
  d = vcat(d1, d2) # the global offset for `self`

  dom_mod = dom[I]
  J = Tuple(i + d) # The index of the codomain of the map to be computed here
  cod_mod = cod[J]
  img_gens = elem_type(cod_mod)[]
  g = gens(dom_mod)
  I_inc = Tuple(-i[1:length(d1)] - d1) # The index of the incoming map of fac.domain_morphism
  I_out = Tuple(i[length(d1)+1:end] + d2) # The index of the outgoing map of fac.codomain_morphism

  if fac.domain_morphism !== nothing && fac.codomain_morphism !== nothing
    img_gens = elem_type(cod_mod)[homomorphism_to_element(cod_mod, compose(compose(fac.domain_morphism[I_inc], element_to_homomorphism(v)), fac.codomain_morphism[I_out])) for v in g]
  elseif fac.domain_morphism !== nothing && fac.codomain_morphism === nothing
    img_gens = elem_type(cod_mod)[homomorphism_to_element(cod_mod, compose(fac.domain_morphism[I_inc], element_to_homomorphism(v))) for v in g]
  elseif fac.domain_morphism === nothing && fac.codomain_morphism !== nothing
    img_gens = elem_type(cod_mod)[homomorphism_to_element(cod_mod, compose(element_to_homomorphism(v), fac.codomain_morphism[I_out])) for v in g]
  end

  return hom(dom_mod, cod_mod, img_gens; check=false) # Set to false eventually
end

function can_compute(fac::HomComplexMorphismFactory, self::AbsHyperComplexMorphism, I::Tuple)
  i = collect(I)
  dom = domain(self)::HomComplex
  dom_dom = domain(dom)::AbsHyperComplex
  dom_cod = codomain(dom)::AbsHyperComplex
  cod = codomain(self)::HomComplex
  cod_dom = domain(cod)::AbsHyperComplex
  cod_cod = codomain(cod)::AbsHyperComplex
  d1 = fac.domain_morphism !== nothing ? offset(fac.domain_morphism) : [0 for i in 1:dim(dom_dom)]
  d2 = fac.codomain_morphism !== nothing ? offset(fac.codomain_morphism) : [0 for i in 1:dim(cod_cod)]
  d = vcat(d1, d2) # the global offset for `self`

  can_compute_index(dom, I) || return false
  J = Tuple(i + d) # The index of the codomain of the map to be computed here
  can_compute_index(cod, J) || return false
  I_inc = Tuple(-i[1:length(d1)] - d1) # The index of the incoming map of fac.domain_morphism
  I_out = Tuple(i[length(d1)+1:end] + d2) # The index of the outgoing map of fac.codomain_morphism
  if fac.domain_morphism !== nothing
    can_compute_index(fac.domain_morphism, I_inc) || return false
  end
  if fac.codomain_morphism !== nothing
    can_compute_index(fac.codomain_morphism, I_out) || return false
  end
  return true
end


@attributes mutable struct HomComplexMorphism{DomainType, CodomainType, MorphismType} <: AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType, HomComplexMorphism{DomainType, CodomainType, MorphismType}}
  internal_morphism::HyperComplexMorphism{DomainType, CodomainType, MorphismType}
  domain_morphism::Union{Nothing, <:AbsHyperComplexMorphism}
  codomain_morphism::Union{Nothing, <:AbsHyperComplexMorphism}

  function HomComplexMorphism(
      ::Type{MorphismType},
      dom::HomComplex,
      cod::HomComplex,
      dom_mor::Union{Nothing, <:AbsHyperComplexMorphism},
      cod_mor::Union{Nothing, <:AbsHyperComplexMorphism}
  ) where {MorphismType}
    dom_dom = domain(dom)::AbsHyperComplex
    dom_cod = codomain(dom)::AbsHyperComplex
    cod_dom = domain(cod)::AbsHyperComplex
    cod_cod = codomain(cod)::AbsHyperComplex

    if dom_mor !== nothing
      @assert dom_dom === codomain(dom_mor)
      @assert cod_dom === domain(dom_mor)
    end

    if cod_mor !== nothing
      @assert dom_cod === domain(cod_mor)
      @assert cod_cod === codomain(cod_mor)
    end

    d1 = dom_mor !== nothing ? offset(dom_mor) : [0 for i in 1:dim(dom_dom)]
    d2 = cod_mor !== nothing ? offset(cod_mor) : [0 for i in 1:dim(cod_cod)]
    d = vcat(d1, d2) # the global offset for this morphism

    map_factory = HomComplexMorphismFactory(MorphismType, dom_mor, cod_mor)
    internal_morphism = HyperComplexMorphism(dom, cod, map_factory, cached=true, offset=d)
    return new{typeof(dom), typeof(cod), MorphismType}(internal_morphism, dom_mor, cod_mor)
  end
end

underlying_morphism(phi::HomComplexMorphism) = phi.internal_morphism

### Specialized functionality
domain_morphism(phi::HomComplexMorphism) = phi.domain_morphism
codomain_morphism(phi::HomComplexMorphism) = phi.codomain_morphism
