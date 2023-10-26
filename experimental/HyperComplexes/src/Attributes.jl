### asking for the entries of the complex
function getindex(HC::HyperComplex, i::Tuple)
  haskey(HC.chains, i) && return HC.chains[i]

  can_compute(HC.chain_factory, HC, i) || error("index out of bounds")

  result = HC.chain_factory(HC, i)
  HC.chains[i] = result
  return result
end

function has_index(HC::HyperComplex, i::Tuple)
  return haskey(HC.chains, i)
end

function can_compute_index(HC::HyperComplex, i::Tuple)
  return can_compute(HC.chain_factory, HC, i)
end

### accessing the maps
function emanating_maps(HC::HyperComplex{ChainType, MorphismType}, i::Tuple) where {ChainType, MorphismType}
  haskey(HC.morphisms, i) && return HC.morphisms[i]
  md = Dict{Int, MorphismType}()
  HC.morphisms = md
  return md
end

function map(HC::HyperComplex, p::Int, i::Tuple)
  phi = emanating_maps(HC, i)
  haskey(phi, p) && return phi[p]
  result = HC.map_factory(HC, p, i)
  phi[p] = result
  return result
end

function has_map(HC::HyperComplex, p::Int, i::Tuple)
  phi = emanating_maps(HC, i)
  return haskey(phi, p)
end

function can_compute_map(HC::HyperComplex, p::Int, i::Tuple)
  return can_compute(HC.map_factory, HC, p, i)
end

### properties
function direction(HC::HyperComplex, p::Int)
  return HC.directions[p]
end

function dim(HC::HyperComplex)
  return HC.d
end

function is_complete(HC::HyperComplex)
  isdefined(HC, :is_complete) && return HC.is_complete
  error("not implemented")
end

function has_upper_bound(HC::HyperComplex, p::Int) 
  return HC.upper_bounds[p] !== nothing
end

function has_lower_bound(HC::HyperComplex, p::Int)
  return HC.lower_bounds[p] !== nothing
end

function upper_bound(HC::HyperComplex, p::Int) 
  return HC.upper_bounds[p]::Int
end

function lower_bound(HC::HyperComplex, p::Int)
  return HC.lower_bounds[p]::Int
end

