  ########################################################################
  # concrete hypercomplexes
  ########################################################################

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
    HC.is_complete === true && return true
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


########################################################################
# concrete double complexes
########################################################################

### Basic getters
function horizontal_direction(D::DoubleComplexOfMorphisms)
  return D.horizontal_direction
end

function vertical_direction(D::DoubleComplexOfMorphisms)
  return D.vertical_direction
end

function chain_cache(D::DoubleComplexOfMorphisms) 
  return D.chains
end

function horizontal_map_cache(D::DoubleComplexOfMorphisms)
  return D.horizontal_maps
end

function vertical_map_cache(D::DoubleComplexOfMorphisms)
  return D.vertical_maps
end

### Boundedness
has_right_bound(D::DoubleComplexOfMorphisms) = D.right_bound !== nothing
has_left_bound(D::DoubleComplexOfMorphisms) = D.left_bound !== nothing
has_upper_bound(D::DoubleComplexOfMorphisms)   = D.upper_bound !== nothing
has_lower_bound(D::DoubleComplexOfMorphisms)   = D.lower_bound !== nothing

function is_complete(D::DoubleComplexOfMorphisms)
  D.is_complete === true && return true

  todo = keys(D.chains)
  isempty(todo) && return false
  for (i, j) in todo
    iszero(D[i, j]) && continue
    has_index(D, i+1, j) || !can_compute_index(D, i+1, j) || return false
    has_index(D, i-1, j) || !can_compute_index(D, i-1, j) || return false
    has_index(D, i, j+1) || !can_compute_index(D, i, j+1) || return false
    has_index(D, i, j-1) || !can_compute_index(D, i, j-1) || return false
  end
  D.is_complete = true
  return true
end

@doc raw"""
    finalize!(D::DoubleComplexOfMorphisms)

Compute all (computable) non-zero entries `D[i, j]` and maps of contiguous pieces of entries 
of `D` in the grid of the double complex, starting from non-zero entries in 
the cache which have already been computed.

This can be used to write a full double complex to disc by filling the cache first.

!!! note This method might not terminate. Use carefully!
"""
function finalize!(D::DoubleComplexOfMorphisms)
  isempty(keys(D.chains)) && error("can not finalize starting from empty cache")
  for (i, j) in keys(D.chains)
    _compute_neighbors_and_maps(D, i, j)
  end
end

function _compute_neighbors_and_maps(D::DoubleComplexOfMorphisms, i::Int, j::Int)
  iszero(D[i, j]) && return
  if !has_index(D, i+1, j) && can_compute_index(D, i+1, j) && !iszero(D[i+1, j])
    _compute_neighbors_and_maps(D, i+1, j) # trigger the computation
  end
  if has_index(D, i+1, j) && !iszero(D[i+1, j])
    if horizontal_direction(D) == :chain && !has_horizontal_map(D, i+1, i) && can_compute_horizontal_map(D, i+1, j)
      horizontal_map(D, i+1, j)
    elseif horizontal_direction(D) == :cochain && !has_horizontal_map(D, i, j) && can_compute_horizontal_map(D, i, j)
      horizontal_map(D, i, j)
    end
  end

  if !has_index(D, i-1, j) && can_compute_index(D, i-1, j) && !iszero(D[i-1, j])
    _compute_neighbors_and_maps(D, i-1, j) # trigger the computation
  end
  if has_index(D, i-1, j) && !iszero(D[i-1, j])
    if horizontal_direction(D) == :cochain && !has_horizontal_map(D, i-1, i) && can_compute_horizontal_map(D, i-1, j)
      horizontal_map(D, i-1, j)
    elseif horizontal_direction(D) == :chain && !has_horizontal_map(D, i, j) && can_compute_horizontal_map(D, i, j)
      horizontal_map(D, i, j)
    end
  end

  if !has_index(D, i, j+1) && can_compute_index(D, i, j+1) && !iszero(D[i, j+1])
    _compute_neighbors_and_maps(D, i, j+1) # trigger the computation
  end
  if has_index(D, i, j+1) && !iszero(D[i, j+1])
    if vertical_direction(D) == :chain && !has_vertical_map(D, i, j+1) && can_compute_vertical_map(D, i, j+1)
      vertical_map(D, i, j+1)
    elseif vertical_direction(D) == :cochain && !has_vertical_map(D, i, j) && can_compute_vertical_map(D, i, j)
      vertical_map(D, i, j)
    end
  end

  if !has_index(D, i, j-1) && can_compute_index(D, i, j-1) && !iszero(D[i, j-1])
    _compute_neighbors_and_maps(D, i, j-1) # trigger the computation
  end
  if has_index(D, i, j-1) && !iszero(D[i, j-1])
    if vertical_direction(D) == :cochain && !has_vertical_map(D, i, j-1) && can_compute_vertical_map(D, i, j-1)
      vertical_map(D, i, j-1)
    elseif vertical_direction(D) == :chain && !has_vertical_map(D, i, j) && can_compute_vertical_map(D, i, j)
      vertical_map(D, i, j)
    end
  end
end


function horizontal_range(D::DoubleComplexOfMorphisms) 
  is_horizontally_bounded(D) || error("complex is not known to be horizontally bounded")
  if horizontal_direction(D) == :chain
    return right_bound(D):-1:left_bound(D)
  elseif horizontal_direction(D) == :cochain
    return left_bound(D):right_bound(D)
  end
  error("typ not recognized")
end

function vertical_range(D::DoubleComplexOfMorphisms) 
  is_vertically_bounded(D) || error("complex is not known to be vertically bounded")
  if vertical_direction(D) == :chain
    return upper_bound(D):-1:lower_bound(D)
  elseif vertical_direction(D) == :cochain
    return lower_bound(D):upper_bound(D)
  end
  error("typ not recognized")
end

function right_bound(D::DoubleComplexOfMorphisms) 
  D.right_bound !== nothing || error("complex has no known horizontal upper bound")
  return D.right_bound::Int
end

function left_bound(D::DoubleComplexOfMorphisms) 
  D.left_bound !== nothing || error("complex has no known horizontal lower bound")
  return D.left_bound::Int
end

function upper_bound(D::DoubleComplexOfMorphisms) 
  D.upper_bound !== nothing || error("complex has no known vertical upper bound")
  return D.upper_bound::Int
end

function lower_bound(D::DoubleComplexOfMorphisms) 
  D.lower_bound !== nothing || error("complex has no known vertical lower bound")
  return D.lower_bound::Int
end

### User facing functionality
function getindex(D::DoubleComplexOfMorphisms, i::Int, j::Int)
  # load from cache if applicable
  haskey(D.chains, (i, j)) && return D.chains[(i, j)]

  # legitimize request
  can_compute_index(D, i, j) || error("index out of bounds")

  # produce entry otherwise
  new_chain = D.chain_factory(D, i, j)
  D.chains[(i, j)] = new_chain

  return new_chain
end

function has_index(D::DoubleComplexOfMorphisms, i::Int, j::Int)
  return haskey(D.chains, (i, j))
end

function can_compute_index(D::DoubleComplexOfMorphisms, i::Int, j::Int)
  has_index(D, i, j) && return true
  return can_compute(D.chain_factory, D, i, j)
end

function vertical_map(D::DoubleComplexOfMorphisms, i::Int, j::Int)
  haskey(vertical_map_cache(D), (i, j)) && return vertical_map_cache(D)[(i, j)]
  can_compute_vertical_map(D, i, j) || error("index out of bounds")
  new_map = D.vertical_map_factory(D, i, j)
  vertical_map_cache(D)[(i, j)] = new_map
  return new_map
end

function has_vertical_map(D::DoubleComplexOfMorphisms, i::Int, j::Int)
  return haskey(D.vertical_maps, (i, j))
end

function can_compute_vertical_map(D::DoubleComplexOfMorphisms, i::Int, j::Int)
  has_vertical_map(D, i, j) && return true
  return can_compute(D.vertical_map_factory, D, i, j)
end

function horizontal_map(D::DoubleComplexOfMorphisms, i::Int, j::Int)
  haskey(horizontal_map_cache(D), (i, j)) && return horizontal_map_cache(D)[(i, j)]
  new_map = D.horizontal_map_factory(D, i, j)
  can_compute_horizontal_map(D, i, j) || error("index out of bounds")
  horizontal_map_cache(D)[(i, j)] = new_map
  return new_map
end

function has_horizontal_map(D::DoubleComplexOfMorphisms, i::Int, j::Int)
  return haskey(D.horizontal_maps, (i, j))
end

function can_compute_horizontal_map(D::DoubleComplexOfMorphisms, i::Int, j::Int)
  has_horizontal_map(D, i, j) && return true
  return can_compute(D.horizontal_map_factory, D, i, j)
end

### Implement the AbsHyperComplex interface for the generic case
getindex(D::DoubleComplexOfMorphisms, i::Tuple) = D[i...]
has_index(D::DoubleComplexOfMorphisms, i::Tuple) = has_index(D, i...)
can_compute_index(D::DoubleComplexOfMorphisms, i::Tuple) = can_compute_index(D, i...)

function map(D::DoubleComplexOfMorphisms, p::Int, i::Tuple)
  isone(p) && return horizontal_map(D, i...)
  p == 2 && return vertical_map(D, i...)
  error("index out of bounds")
end

function has_map(D::DoubleComplexOfMorphisms, p::Int, i::Tuple)
  isone(p) && return has_horizontal_map(D, i...)
  p == 2 && return has_vertical_map(D, i...)
  error("index out of bounds")
end

function can_compute_map(D::DoubleComplexOfMorphisms, p::Int, i::Tuple)
  isone(p) && return can_compute_horizontal_map(D, i...)
  p == 2 && return can_compute_vertical_map(D, i...)
  error("index out of bounds")
end

function direction(D::DoubleComplexOfMorphisms, p::Int)
  isone(p) && return horizontal_direction(D)
  p == 2 && return vertical_direction(D)
  error("index out of bounds")
end

dim(D::DoubleComplexOfMorphisms) = 2

function has_upper_bound(D::DoubleComplexOfMorphisms, p::Int)
  isone(p) && return has_right_bound(D)
  p == 2 && return has_upper_bound(D)
  error("index out of bounds")
end

function has_lower_bound(D::DoubleComplexOfMorphisms, p::Int)
  isone(p) && return has_left_bound(D)
  p == 2 && return has_lower_bound(D)
  error("index out of bounds")
end

function upper_bound(D::DoubleComplexOfMorphisms, p::Int)
  isone(p) && return right_bound(D)
  p == 2 && return upper_bound(D)
  error("index out of bounds")
end

function lower_bound(D::DoubleComplexOfMorphisms, p::Int)
  isone(p) && return left_bound(D)
  p == 2 && return lower_bound(D)
  error("index out of bounds")
end

