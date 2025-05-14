struct ResolutionModuleFactory{ChainType, MapType} <: HyperComplexChainFactory{ChainType}
  orig_mod::SubquoModule
  map_cache::Vector{MapType}

  function ResolutionModuleFactory(M::SubquoModule)
    R = base_ring(M)
    ChainType = FreeMod{elem_type(R)}
    MapType = FreeModuleHom{ChainType, ChainType, Nothing}
    map_cache = MapType[]
    return new{ChainType, MapType}(M, map_cache)
  end
end

# Generic Code via iterated kernel computation
function (fac::ResolutionModuleFactory{ChainType})(c::AbsHyperComplex, I::Tuple) where {ChainType <: SparseFPModule}
  i = first(I)
  R = base_ring(fac.orig_mod)

  if iszero(i)
    return _make_free_module(fac.orig_mod, gens(fac.orig_mod))
  end

  if isone(i)
    aug = hom(c[0], fac.orig_mod, gens(fac.orig_mod); check=false)
    aug.generators_map_to_generators = true
    K, inc = kernel(aug)
    next = _make_free_module(K, gens(K))
    phi = hom(next, c[0], ambient_representatives_generators(K); check=false)
    push!(fac.map_cache, phi)
    return next
  end

  prev = c[i-1]
  if iszero(prev)
    next, inc = zero_object(prev)
    push!(fac.map_cache, inc)
    return next
  end

  phi = map(c, 1, (i-1,))
  K, inc = kernel(phi)
  if iszero(K)
    next, inc = zero_object(prev)
    push!(fac.map_cache, inc)
    return next
  end
  next = _make_free_module(K, gens(K))
  phi = hom(next, c[i-1], ambient_representatives_generators(K); check=false)
  push!(fac.map_cache, phi)
  
  return next
end

# Modified version to use Schreyer's resolution where applicable
function (fac::ResolutionModuleFactory{ChainType})(
             c::AbsHyperComplex, I::Tuple
           ) where {ChainType <: SparseFPModule{<:MPolyRingElem{<:FieldElem}}}
  i = first(I)
  R = base_ring(fac.orig_mod)

  is_zero(i) && return _make_free_module(fac.orig_mod, gens(fac.orig_mod))

  if is_empty(fac.map_cache)
    # start the resolution from scratch
    # It turns out to usually be better to compute the full Schreyer resolution at once 
    # instead of incrementally.
    res = free_resolution(fac.orig_mod; algorithm=:fres)
    phi = map(res, 1)
    F = c[0] # caught above
    id = hom(codomain(phi), F, gens(F))
    push!(fac.map_cache, compose(phi, id))
    r = first(range(res.C)) # the index of the last module
    for j in 2:r
      push!(fac.map_cache, map(res, j))
    end
    @assert is_zero(domain(last(fac.map_cache)))
  end

  m = length(fac.map_cache)
  i == 0 && return codomain(first(fac.map_cache))
  i <= m && return domain(fac.map_cache[i])

  for j in m+1:i
    next, inc = zero_object(domain(last(fac.map_cache)))
    push!(fac.map_cache, inc)
  end
  return domain(last(fac.map_cache))
end


function zero_object(M::SparseFPModule)
  if is_graded(M)
    result = graded_free_module(base_ring(M), [])
    return result, hom(result, M, elem_type(M)[]; check=false)
  else
    result = FreeMod(base_ring(M), 0)
    return result, hom(result, M, elem_type(M)[]; check=false)
  end
end


# TODO: Faster code using Schreyer's algorithm for the special case of polynomial rings.

function can_compute(fac::ResolutionModuleFactory, c::AbsHyperComplex, I::Tuple)
  return first(I) >= 0
end

struct ResolutionMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType} end

function can_compute(fac::ResolutionMapFactory, c::AbsHyperComplex, p::Int, I::Tuple)
  isone(p) || return false
  i = first(I)
  return i > 0
end

function (fac::ResolutionMapFactory)(c::AbsHyperComplex, p::Int, I::Tuple)
  chain_fac = chain_factory(c)
  i = first(I)
  c[i] # Fill the cache
  return chain_fac.map_cache[i]
end


function free_resolution(::Type{T}, M::SubquoModule{RET}) where {T<:SimpleFreeResolution, RET}
  ChainType = FreeMod{RET}
  MorphismType = FreeModuleHom{ChainType, ChainType, Nothing}

  chain_fac = ResolutionModuleFactory(M)
  map_fac = ResolutionMapFactory{MorphismType}()

  R = base_ring(M)
  upper_bound = (R isa MPolyRing ? ngens(R) : nothing)
  internal_complex = HyperComplex(1, chain_fac, map_fac, [:chain],
                                  upper_bounds = [upper_bound], 
                                  lower_bounds = [0]
                                 )
  result = SimpleFreeResolution(M, internal_complex)
  MC = ZeroDimensionalComplex(M)[0:0] # Wrap MC as a 1-dimensional complex concentrated in degree 0
  aug_map = hom(result[(0,)], M, gens(M); check=false) # The actual augmentation map
  aug_map.generators_map_to_generators = true
  aug_map_comp = MorphismFromDict(result, MC, Dict{Tuple, typeof(aug_map)}([(0,)=>aug_map]))
  result.augmentation_map = aug_map_comp
  return result, aug_map_comp
end

function free_resolution(::Type{T}, I::Ideal{RET}) where {T<:SimpleFreeResolution, RET}
  R = base_ring(I)
  F = (!is_graded(R) ? FreeMod(R, 1) : graded_free_module(R, [zero(grading_group(R))]))
  M, _ = I*F
  if is_graded(R)
    @assert is_graded(M)
  end
  return free_resolution(T, M)
end

function free_resolution(::Type{T}, A::MPolyQuoRing) where {T<:SimpleFreeResolution}
  R = base_ring(A)
  I = modulus(A)
  F = (!is_graded(R) ? FreeMod(R, 1) : graded_free_module(R, [zero(grading_group(R))]))
  M, _ = I*F
  AA, _ = quo(F, M)
  if is_graded(R)
    @assert is_graded(M)
  end
  return free_resolution(T, AA)
end

### Additional getters
augmentation_map(c::SimpleFreeResolution) = c.augmentation_map

### Additional functionality
function betti(b::SimpleFreeResolution; project::Union{FinGenAbGroupElem, Nothing} = nothing, reverse_direction::Bool = false)
  return betti_table(b; project, reverse_direction)
end
function betti_table(C::AbsHyperComplex; 
    lower_bound::Union{Int, Nothing}=(has_lower_bound(C) ? Oscar.lower_bound(C) : nothing),
    upper_bound::Union{Int, Nothing}=(has_upper_bound(C) ? Oscar.upper_bound(C) : nothing),
    project::Union{FinGenAbGroupElem, Nothing} = nothing, reverse_direction::Bool=false
  )
  @assert dim(C) == 1 "complex must be one-dimensional"
  @assert lower_bound !== nothing && upper_bound !== nothing "explicit bounds must be known"
  @assert direction(C, 1) == :chain "only implemented for chain complexes"
  generator_count = Dict{Tuple{Int, Any}, Int}()
  for i in lower_bound:upper_bound
    @assert is_graded(C[i]) "one of the modules in the graded free resolution is not graded"
    module_degrees = degree.(gens(C[i]))
    for degree in module_degrees
      idx = (i, degree)
      generator_count[idx] = get(generator_count, idx, 0) + 1
    end
  end
  return BettiTable(generator_count, project = project, reverse_direction = reverse_direction)
end

function minimal_betti_table(C::AbsHyperComplex;
    lower_bound::Union{Int, Nothing}=(has_lower_bound(C) ? Oscar.lower_bound(C) : nothing),
    upper_bound::Union{Int, Nothing}=(has_upper_bound(C) ? Oscar.upper_bound(C) : nothing)
  )
  @assert dim(C) == 1 "complex must be one-dimensional"
  @assert lower_bound !== nothing && upper_bound !== nothing "explicit bounds must be known"
  @assert direction(C, 1) == :chain "only implemented for chain complexes"
  offsets = Dict{FinGenAbGroupElem, Int}()
  betti_hash_table = Dict{Tuple{Int, Any}, Int}()
  for i in lower_bound+1:upper_bound+1
    phi = map(C, i)
    F = domain(phi)
    @assert is_graded(F) "modules must be graded"
    @assert is_standard_graded(base_ring(F)) "ring must be standard graded"
    G = codomain(phi)
    @assert is_graded(G) "modules must be graded"
    @assert is_standard_graded(base_ring(G)) "ring must be standard graded"
    dom_degs = unique!([degree(g) for g in gens(F)])
    cod_degs = unique!([degree(g) for g in gens(G)])
    for d in cod_degs
      d::FinGenAbGroupElem
      if d in dom_degs
        _, _, sub_mat = _constant_sub_matrix(phi, d)
        r = rank(sub_mat)
        c = ncols(sub_mat) - r - get(offsets, d, 0)
        !iszero(c) && (betti_hash_table[(i-1, d)] = c)
        offsets[d] = r
      else
        c = length(_indices_of_generators_of_degree(G, d)) - get(offsets, d, 0)
        !iszero(c) && (betti_hash_table[(i-1, d)] = c)
      end
    end
  end
  return BettiTable(betti_hash_table)
end
