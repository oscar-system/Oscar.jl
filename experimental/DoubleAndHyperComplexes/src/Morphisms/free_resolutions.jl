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
function (fac::ResolutionModuleFactory{ChainType})(c::AbsHyperComplex, I::Tuple) where {ChainType <: ModuleFP}
  i = first(I)
  R = base_ring(fac.orig_mod)

  if iszero(i)
    n = ngens(fac.orig_mod)
    return _make_free_module(fac.orig_mod, gens(fac.orig_mod))
  end

  if isone(i)
    aug = hom(c[0], fac.orig_mod, gens(fac.orig_mod))
    K, inc = kernel(aug)
    next = _make_free_module(K, gens(K))
    phi = hom(next, c[0], ambient_representatives_generators(K))
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
  phi = hom(next, c[i-1], ambient_representatives_generators(K))
  push!(fac.map_cache, phi)
  
  return next
end

function zero_object(M::ModuleFP)
  if is_graded(M)
    result = graded_free_module(base_ring(M), [])
    return result, hom(result, M, elem_type(M)[])
  else
    result = FreeMod(base_ring(M), 0)
    return result, hom(result, M, elem_type(M)[])
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
  aug_map = hom(result[(0,)], M, gens(M)) # The actual augmentation map
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

### Additional getters
augmentation_map(c::SimpleFreeResolution) = c.augmentation_map
