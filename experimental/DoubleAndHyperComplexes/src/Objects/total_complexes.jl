### Production of the chains in a total complex
struct TotalComplexChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  orig::AbsHyperComplex

  index_cache::Dict{Int, Vector{<:Tuple}} # index_cache[k] contains all tuples of indices 
                                        # with non-trivial contributions to the direct sum
                                        # in degree k.
  injections::Dict{Int, Vector{<:Map}}
  projections::Dict{Int, Vector{<:Map}}
  function TotalComplexChainFactory(c::AbsHyperComplex{ChainType}) where {ChainType}
    return new{ChainType}(c, Dict{Int, Tuple}(), Dict{Int, Vector{Map}}(), Dict{Int, Vector{Map}}())
  end
end

original_complex(fac::TotalComplexChainFactory) = fac.orig
index_cache(fac::TotalComplexChainFactory) = fac.index_cache
injections(fac::TotalComplexChainFactory) = fac.injections
projections(fac::TotalComplexChainFactory) = fac.projections

function (fac::TotalComplexChainFactory{ChainType})(c::AbsHyperComplex, I::Tuple) where {ChainType}
  d = first(I)
  new_indices = Vector{Tuple}()
  summands = ChainType[]
  orig = original_complex(fac)

  # determine a range over which to iterate
  if all(k->has_lower_bound(orig, k), 1:dim(orig))
    b0 = [lower_bound(orig, k) for k in 1:dim(orig)]
    b = sum(b0; init=0)
    for j in weak_compositions(d - b, dim(orig))
      j = j + b0
      any(k->(has_upper_bound(orig, k) && j[k] > upper_bound(orig, k)), 1:dim(orig)) && continue
      J = Tuple(j)
      if can_compute_index(orig, J)
        M = orig[J]
        push!(summands, M)
        push!(new_indices, J)
      end
    end
  else # all(k->has_upper_bound(orig, k), 1:dim(orig)) must be true then
    b0 = [upper_bound(orig, k) for k in 1:dim(orig)]
    b = sum(b0; init=0)
    for j in weak_compositions(b - d, dim(orig))
      j = b0 - j
      any(k->(has_lower_bound(orig, k) && j[k] < lower_bound(orig, k)), 1:dim(orig)) && continue
      J = Tuple(j)
      if can_compute_index(orig, J)
        M = orig[J]
        push!(summands, M)
        push!(new_indices, J)
      end
    end
  end
  isempty(summands) && error("can not compute the direct some of an empty list")
  tmp_res = _direct_sum(summands)
  @assert tmp_res isa Tuple{<:ChainType, <:Vector, <:Vector} "the output of `direct_sum` does not have the anticipated format; see the source code for details"
  # If you got here because of the error message above:
  # This is supposed to be generic code and it attempts to build the direct sum 
  # of `summands` using the method of `direct_sum` for your specific type. 
  # The convention in Oscar is that this should return a triple `(S, inj, pr)` consisting 
  # of the actual direct sum `S`, a `Vector` of injections `inj` and a `Vector` of 
  # projections `pr` 
  # Unfortunately, it can not be assumed that the original method for `direct_sum` 
  # produces this output and it is an ongoing effort to streamline this throughout 
  # Oscar. Since such changes tend to happen with severe delay, if ever, we provide 
  # a little workaround here. If this code does not run for your type of chains, 
  # you may try two things:
  #
  #   1) Overwrite the method for `_direct_sum` below for your type and wrap the 
  #      original method for `direct_sum` so that the anticipated output is produced
  #
  #   2) If that does not work, you may also overwrite the whole method for production 
  #      of the direct sum here.
  result, inj, pr = tmp_res
  index_cache(fac)[d] = new_indices
  injections(fac)[d] = Map[i for i in inj]
  projections(fac)[d] = Map[x for x in pr]
  return result
end

@doc raw"""
    _direct_sum(u::Vector{T}) where {T}

Internal method to return a triple `(s, incs, prs)` consisting of an object 
`s` representing the direct sum of the entries in `u`, together with vectors 
of maps `incs` for the inclusion maps `u[i] → s` and `prs` for 
the projections `s →  u[i]`.

Generically this will default to `direct_sum(u)`. If that does not produce 
a result with the required output format, you must overwrite this method 
for your specific type `T`.
"""
function _direct_sum(u::Vector{T}) where {T}
  return direct_sum(u...)
end

# overwriting the method for finitely generated modules
function _direct_sum(u::Vector{T}) where {T<:ModuleFP}
  return direct_sum(u...; task=:both)
end


### Production of the maps of a total complex
struct TotalComplexMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType}
  orig::AbsHyperComplex

  function TotalComplexMapFactory(c::AbsHyperComplex{<:Any, MorphismType}) where {MorphismType} 
    return new{MorphismType}(c)
  end
end

original_complex(fac::TotalComplexMapFactory) = fac.orig

function (fac::TotalComplexMapFactory)(c::AbsHyperComplex, p::Int, I::Tuple)
  d = first(I)
  orig = original_complex(fac)
  chain_fac = chain_factory(c)
  inc = (direction(c, 1) == :chain ? -1 : 1)
  next = d + inc
  dom = c[d]
  cod = c[next]
  result = hom(dom, cod, elem_type(cod)[zero(cod) for i in 1:ngens(dom)]; check=false)
  for (ind, J) in enumerate(index_cache(chain_fac)[d])
    for k in 1:dim(orig)
      target = collect(J) + (direction(orig, k) == :chain ? -1 : 1)*[(l == k ? 1 : 0) for l in 1:dim(orig)]
      T = Tuple(target)
      index_in_cod = findfirst(==(T), index_cache(chain_fac)[next])
      index_in_cod === nothing && continue
      phi = map(orig, k, J)
      @assert codomain(phi) === orig[T]
      inc_map = compose(compose(projections(chain_fac)[d][ind], phi), 
                        injections(chain_fac)[next][index_in_cod])
      # We want the following sign convention for the total complex 
      # of a double complex:
      #
      #       •  ←   •  ←  •  ←  •  ←  •  ←  •  ←
      #       ↓    -1↓     ↓   -1↓     ↓   -1↓      
      #       •  ←   •  ←  •  ←  •  ←  •  ←  •  ←
      #       ↓    -1↓     ↓   -1↓     ↓   -1↓      
      #       •  ←   •  ←  •  ←  •  ←  •  ←  •  ←
      #       ↓    -1↓     ↓   -1↓     ↓   -1↓      
      #       •  ←   •  ←  •  ←  •  ←  •  ←  •  ←
      #
      # Inductively we get a formula for the sign convention in 
      # the total complex of an arbitrary hyper-complex which is 
      # implemented below. To this end, we consider tot(c) as 
      # tot(d) where d is the double complex whose columns consist of the 
      # total complexes of the slices of c for fixed indices in its first 
      # dimension and with the horizontal maps of d the induced maps 
      # in the total complexes of the slices.
      result = result + (-1)^(sum(J[k+1:end]; init=0)) * inc_map
    end
  end
  return result
end

function can_compute(fac::TotalComplexMapFactory, c::AbsHyperComplex, p::Int, i::Tuple)
  I = collect(i)
  next = I + (direction(c, p) == :chain ? -1 : 1)*[k==p ? 1 : 0 for k in 1:dim(c)]
  return can_compute_index(c, i) && can_compute_index(c, Tuple(next))
end

function can_compute(fac::TotalComplexChainFactory{ChainType}, c::AbsHyperComplex, I::Tuple) where {ChainType}
  # We must find at least one computable entry in this slice
  d = first(I)
  new_indices = Vector{Tuple}()
  summands = ChainType[]
  orig = original_complex(fac)

  # determine a range over which to iterate
  if all(k->has_lower_bound(orig, k), 1:dim(orig))
    b0 = [lower_bound(orig, k) for k in 1:dim(orig)]
    b = sum(b0; init=0)
    d - b < 0 && return false
    for j in weak_compositions(d - b, dim(orig))
      j = j + b0
      any(k->(has_upper_bound(orig, k) && j[k] > upper_bound(orig, k)), 1:dim(orig)) && continue
      J = Tuple(j)
      can_compute_index(orig, J) && return true
    end
  else # all(k->has_upper_bound(orig, k), 1:dim(orig)) must be true then
    b0 = [upper_bound(orig, k) for k in 1:dim(orig)]
    b = sum(b0; init=0)
    b - d < 0 && return false
    for j in weak_compositions(b - d, dim(orig))
      j = b0 - j
      J = Tuple(j)
      can_compute_index(orig, J) && return true
    end
  end
  return false
end

function total_complex(c::AbsHyperComplex)
  upper = (all(k->has_upper_bound(c, k), 1:dim(c)) ? sum(upper_bound(c, k) for k in 1:dim(c); init=0) : nothing)
  lower = (all(k->has_lower_bound(c, k), 1:dim(c)) ? sum(lower_bound(c, k) for k in 1:dim(c); init=0) : nothing)
  all(k->direction(c, k) == :chain, 1:dim(c)) || all(k->direction(c, k) == :cochain, 1:dim(c)) || error("complex must be all chain or all cochain")
  dir = (all(k->direction(c, k) == :chain, 1:dim(c)) ? :chain : :cochain)
  chain_fac = TotalComplexChainFactory(c)
  map_fac = TotalComplexMapFactory(c)
  complex = HyperComplex(1, chain_fac, map_fac, [dir],
                         upper_bounds = [upper], 
                         lower_bounds = [lower]
                        )
  return TotalComplex(c, complex)
end

original_complex(tot::TotalComplex) = tot.original

# Returns the list of indices l in the `original_complex` 
# of `tot` which make up the degree i direct sum.
function indices_in_summand(tot::TotalComplex, i::Int)
  has_index(tot, (i,))|| tot[i] # Fill the cache
  return index_cache(chain_factory(tot))[i]
end

function injections_for_summand(tot::TotalComplex, i::Int)
  has_index(tot, (i,))|| tot[i] # Fill the cache
  return injections(chain_factory(tot))[i]
end

function projections_for_summand(tot::TotalComplex, i::Int)
  has_index(tot, (i,))|| tot[i] # Fill the cache
  return projections(chain_factory(tot))[i]
end

# TODO: Cache these and use dictionaries
function injection(tot::TotalComplex, i::Tuple)
  d = sum(i)
  v = indices_in_summand(tot, d)
  k = findfirst(==(i), v)
  k === nothing && return nothing
  return injections_for_summand(tot, d)[k]
end

function projection(tot::TotalComplex, i::Tuple)
  d = sum(i)
  v = indices_in_summand(tot, d)
  k = findfirst(==(i), v)
  k === nothing && return nothing
  return projections_for_summand(tot, d)[k]
end

homology(tot::TotalComplex, i::Int) = homology(underlying_complex(tot), 1, (i,))
kernel(tot::TotalComplex, i::Int) = kernel(underlying_complex(tot), 1, (i,))
boundary(tot::TotalComplex, i::Int) = boundary(underlying_complex(tot), 1, (i,))

