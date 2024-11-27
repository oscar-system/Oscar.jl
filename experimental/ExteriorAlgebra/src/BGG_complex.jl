# Helper structure to transfer back and forth between S- and E-modules
mutable struct BGGHelper{T}
  E::ExteriorAlgebra{T}
  S::MPolyDecRing{T}
  d::Int

  # Fields for caching
  ind_to_exp_dom::Vector{Vector{Int}}
  exp_to_ind_dom::Dict{Vector{Int}, Int}
  ind_to_exp_cod::Vector{Vector{Int}}
  exp_to_ind_cod::Dict{Vector{Int}, Int}

  dom::FreeMod{ExtAlgElem{T}} # the (co)domains of the BGG morphism
  cod::FreeMod{ExtAlgElem{T}}
  phi::FreeModuleHom{FreeMod{ExtAlgElem{T}}, FreeMod{ExtAlgElem{T}}, Nothing} # the actual morphism
  K::SubquoModule{ExtAlgElem{T}} # its kernel
  inc::SubQuoHom # the inclusion K ↪ dom
  res::SimpleFreeResolution # a free resolution P of K
  aug::FreeModuleHom # the augmentation map P₀ →  K
  all_monomials::Vector{<:MPolyDecRingElem}

  function BGGHelper(E::ExteriorAlgebra{T}, S::MPolyDecRing{T}, d::Int) where {T}
    return new{T}(E, S, d)
  end
end

function domain(h::BGGHelper)
  if !isdefined(h, :dom)
    if h.d < 0
      h.dom = graded_free_module(h.E, 0)
    else
      n = rank(h.E)
      h.dom = graded_free_module(h.E, [-h.d for i in 1:length(WeakCompositions(h.d, n))])
    end
  end
  return h.dom
end

function codomain(h::BGGHelper)
  if !isdefined(h, :cod)
    if h.d < -1
      h.cod = graded_free_module(h.E, 0)
    else
      n = rank(h.E)
      h.cod = graded_free_module(h.E, [-(h.d+1) for i in 1:length(WeakCompositions(h.d + 1, n))])
    end
  end
  return h.cod
end

function ind_to_exp_dom(h::BGGHelper)
  if !isdefined(h, :ind_to_exp_dom)
    h.ind_to_exp_dom = [a.c for a in WeakCompositions(h.d, rank(h.E))]
  end
  return h.ind_to_exp_dom
end

function exp_to_ind_dom(h::BGGHelper)
  if !isdefined(h, :exp_to_ind_hom)
    iter = WeakCompositions(h.d, rank(h.E))
    r = length(iter)
    h.exp_to_ind_dom = Dict{Vector{Int}, Int}(a.c => i for (i, a) in zip(1:r, iter))
  end
  return h.exp_to_ind_dom
end

function ind_to_exp_cod(h::BGGHelper)
  if !isdefined(h, :ind_to_exp_cod)
    h.ind_to_exp_cod = [a.c for a in WeakCompositions(h.d+1, rank(h.E))]
  end
  return h.ind_to_exp_cod
end

function exp_to_ind_cod(h::BGGHelper)
  if !isdefined(h, :exp_to_ind_cod)
    iter = WeakCompositions(h.d+1, rank(h.E))
    r = length(iter)
    h.exp_to_ind_cod = Dict{Vector{Int}, Int}(a.c => i for (i, a) in zip(1:r, iter))
  end
  return h.exp_to_ind_cod
end

# take a homogeneous polynomial in S and return the SRow corresponding 
# element in the domain/codomain of the E-modules
function (h::BGGHelper{T})(f::MPolyDecRingElem{T}) where {T}
  degf = Int(degree(f; check=false)[1])
  E = h.E
  S = h.S
  R = base_ring(S)
  if degf == h.d # element of the domain
    to_ind_dom = exp_to_ind_dom(h)
    return sparse_row(R, [to_ind_dom[e] for e in exponents(f)],
                      collect(coefficients(f)))
  elseif degf == h.d + 1
    to_ind_cod = exp_to_ind_cod(h)
    return sparse_row(R, [to_ind_cod[e] for e in exponents(f)],
                      collect(coefficients(f)))
  end
  error("execution should never get here")
end

# take an element in the domain/codomain of h and turn it into a homogeneous polynomial
function (h::BGGHelper)(v::FreeModElem)
  E = h.E
  S = h.S
  n = rank(E)
  if parent(v) === domain(h)
    return h(coordinates(v), :domain)
  elseif parent(v) === codomain(h)
    return h(coordinates(v), :codomain)
  end
  error("execution should never get here")
end

# for a list of integers indicating basis vectors, 
# produce the list of corresponding monomials
function all_monomials(h::BGGHelper)
  S = h.S
  h.d < 0 && return elem_type(S)[]
  R = base_ring(S)
  if !isdefined(h, :all_monomials)
    exp_list = ind_to_exp_dom(h)
    result = sizehint!(elem_type(h.S)[], sizeof(exp_list))
    for e in exp_list
      ctx = MPolyBuildCtx(h.S)
      push_term!(ctx, one(R), e)
      push!(result, finish(ctx))
    end
    h.all_monomials = result
  end
  return h.all_monomials
end

function (h::BGGHelper)(v::SRow, s::Symbol=:domain)
  @assert base_ring(v) === base_ring(h.S)
  ctx = MPolyBuildCtx(h.S)
  if s == :domain
    to_exp_dom = ind_to_exp_dom(h)
    for (i, c) in v
      push_term!(ctx, c, to_exp_dom[i])
    end
    return finish(ctx)
  elseif s == :codomain
    to_exp_cod = ind_to_exp_cod(h)
    for (i, c) in v
      push_term!(ctx, c, to_exp_cod[i])
    end
    return finish(ctx)
  end
  error("execution should never get here")
end

function morphism(h::BGGHelper)
  G = codomain(h)
  img_gens = elem_type(G)[]
  n = rank(h.E)
  h.d < 0 && return hom(domain(h), codomain(h), elem_type(codomain(h))[])
  # iterate through the monomials of degree d in n variables
  for (i, a) in enumerate(WeakCompositions(h.d, n))
    e = a.c # the exponent vector of the monomial
    y = [copy(e) for i in 1:n]
    for j in 1:n
      y[j][j] += 1 # multiplication by xⱼ
    end
    to_ind_cod = exp_to_ind_cod(h)
    indices = [to_ind_cod[a] for a in y]
    push!(img_gens, G(sparse_row(h.E, indices, gens(h.E))))
  end
  return hom(domain(h), G, img_gens)
end

function kernel(h::BGGHelper)
  if !isdefined(h, :K)
    K, inc = kernel(morphism(h))
    h.K = K
    h.inc = inc
  end
  return h.K
end

function kernel_inclusion(h::BGGHelper)
  if !isdefined(h, :inc)
    kernel(h)
  end
  return h.inc
end

function resolution(h::BGGHelper)
  if !isdefined(h, :res)
    res, aug = free_resolution(SimpleFreeResolution, kernel(h))
    h.res = res
    h.aug = aug[0]
  end
  return h.res
end

function augmentation_map(h::BGGHelper)
  if !isdefined(h, :aug)
    resolution(h)
  end
  return h.aug
end
  

### Production of the chains
struct BGGChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  # Fields needed for production
  orig::AbsHyperComplex
  S::MPolyDecRing
  E::ExteriorAlgebra
  d::Int
  helper::Dict{Int, <:BGGHelper}
  helper_combinations::Dict{<:Tuple, <:Vector{<:BGGHelper}}
  bgg_maps::Dict{Tuple, FreeModuleHom}

  function BGGChainFactory(
      S::MPolyDecRing, E::ExteriorAlgebra, orig::AbsHyperComplex{T}, d::Int
    ) where {T<:ModuleFP{<:MPolyDecRingElem}}
    @assert dim(orig) == 1 "complex must be one-dimensional"
    @assert direction(orig, 1) == :chain "only chain complexes are supported"
    @assert is_z_graded(S) "ring must be standard graded"
    return new{ModuleFP{elem_type(E)}}(orig, S, E, d, 
                                       Dict{Int, BGGHelper}(),
                                       Dict{Tuple, Vector{<:BGGHelper}}(),
                                       Dict{Tuple, FreeModuleHom}()
                                      )
  end
end

function get_helper(fac::BGGChainFactory, d::Int)
  if !haskey(fac.helper, d)
    fac.helper[d] = BGGHelper(fac.E, fac.S, d)
  end
  return fac.helper[d]
end

function (fac::BGGChainFactory)(self::AbsHyperComplex, I::Tuple)
  return _build_BGG_module!(fac, fac.orig[I], I)
end

function _build_BGG_module!(fac::BGGChainFactory, F::SubquoModule, I::Tuple)
  error("not implemented")
end

function _build_BGG_module!(fac::BGGChainFactory, F::FreeMod, I::Tuple)
  E = fac.E
  is_zero(rank(F)) && return graded_free_module(E, 0)
  S = fac.S
  n = rank(E)
  s = fac.d
  R = base_ring(E)
  @assert ngens(S) == rank(E)
  @assert base_ring(S) === R
  @assert is_graded(F)
  maps = FreeModuleHom[]
  kernels = SubquoModule[]
  helpers = BGGHelper[]
  for (i, d) in enumerate(degrees_of_generators(F))
    dd = Int(d[1])
    helper = get_helper(fac, s-dd)
    phi = morphism(helper)
    dom_twist, amb_twist = twist(domain(phi), dd)
    phi_twist = twist(phi, dd; domain_twist=amb_twist)
    K_twist, _ = twist(kernel(helper), d; ambient_twist=amb_twist)
    @assert all(Int(e[1]) == -s for e in degrees_of_generators(dom_twist))
    push!(maps, phi_twist)
    push!(helpers, helper)
    push!(kernels, K_twist)
  end
  phi = direct_sum(maps)
  dom = domain(phi)
  cod = codomain(phi)
  fac.helper_combinations[I] = helpers
  fac.bgg_maps[I] = phi
  K = direct_sum(kernels...)[1]
  return K
end

function can_compute(fac::BGGChainFactory, self::AbsHyperComplex, i::Tuple)
  return can_compute_index(fac.orig, i)
end

### Production of the morphisms 
struct BGGMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType} end

function (fac::BGGMapFactory)(self::AbsHyperComplex, p::Int, I::Tuple)
  @assert isone(p)
  i = first(I)
  cfac = chain_factory(self)
  orig = cfac.orig
  dom = self[I]
  cod = self[(i-1,)]
  (is_zero(dom) || is_zero(cod)) && return zero_morphism(dom, cod)
  dom_amb = ambient_free_module(dom)
  cod_amb = ambient_free_module(cod)
  orig_dom = orig[I]::FreeMod
  orig_cod = orig[(i-1,)]::FreeMod
  orig_phi = map(orig, 1, I)
  
  # We first build the morphism on the ambient modules

  # transfer these elements to homogeneous S-module elements
  helpers = cfac.helper_combinations[I]::Vector{<:BGGHelper}
  img_gens = elem_type(orig_dom)[]
  w = zero(orig_dom)
  for (j, h) in enumerate(helpers) # corresponding to the direct summands
    delta = rank(domain(h))
    g = orig_dom[j]
    # the generators of the E-module for this summand were just the 
    # monomials of the S-module in this degree
    img_gens = vcat(img_gens, [x*g for x in all_monomials(h)])
  end

  # now we map the monomials over to the other side
  img_gens2 = orig_phi.(img_gens) 

  # and translate them back
  helpers = cfac.helper_combinations[(i-1,)]
  img_gens3 = elem_type(cod_amb)[]
  offsets = rank.(domain.(helpers))
  for j = length(offsets):-1:2
    offsets[j] = sum(offsets[1:j])
  end
  pushfirst!(offsets, 0)

  for v in img_gens2
    w = zero(cod_amb)
    for (j, poly) in coordinates(v)
      helper = helpers[j]
      coord = helper(poly) # convert the polynomial into an element of the summand
      coord.pos.+=offsets[j] # a dirty hack to shift the vector to its correct position in the total sum
      w += cod_amb(map_entries(cfac.E, coord))
    end
    push!(img_gens3, w)
  end

  # this gives us the map on the ambient free E-modules
  amb_map = hom(dom_amb, cod_amb, img_gens3)

  # finally we need to compute the induced map on the kernels

  img_gens4 = elem_type(cod)[SubquoModuleElem(amb_map(g), cod) for g in ambient_representatives_generators(dom)]
  return hom(dom, cod, img_gens4; check=false)
end

function can_compute(fac::BGGMapFactory, self::AbsHyperComplex, p::Int, i::Tuple)
  return can_compute_map(chain_factory(self).orig, p, i)
end

### The concrete struct
@attributes mutable struct BGGComplex{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType} 
  internal_complex::HyperComplex{ChainType, MorphismType}

  function BGGComplex(E::ExteriorAlgebra{T}, S::MPolyDecRing{T}, orig::AbsHyperComplex{CT, MT}, s::Int) where {T <: RingElem, CT <: ModuleFP, MT<:ModuleFPHom}
    chain_fac = BGGChainFactory(S, E, orig, s)
    map_fac = BGGMapFactory{ModuleFPHom}()

    # Assuming d is the dimension of the new complex
    internal_complex = HyperComplex(1, chain_fac, map_fac, [:chain], 
                                    lower_bounds = Union{Int, Nothing}[has_lower_bound(orig, 1) ? lower_bound(orig, 1) : nothing],
                                    upper_bounds = Union{Int, Nothing}[has_upper_bound(orig, 1) ? upper_bound(orig, 1) : nothing])
    # Assuming that ChainType and MorphismType are provided by the input
    return new{ModuleFP{ExtAlgElem{T}}, ModuleFPHom}(internal_complex)
  end
end

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::BGGComplex) = c.internal_complex

