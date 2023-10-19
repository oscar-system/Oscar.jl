### Basic getters
function horizontal_typ(D::DoubleComplexOfMorphisms)
  return D.horizontal_typ
end

function vertical_typ(D::DoubleComplexOfMorphisms)
  return D.vertical_typ
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
has_right_bound(D::DoubleComplexOfMorphisms) = isdefined(D, :right_bound)
has_left_bound(D::DoubleComplexOfMorphisms) = isdefined(D, :left_bound)
has_upper_bound(D::DoubleComplexOfMorphisms)   = isdefined(D, :upper_bound)
has_lower_bound(D::DoubleComplexOfMorphisms)   = isdefined(D, :lower_bound)

function horizontal_range(D::DoubleComplexOfMorphisms) 
  is_horizontally_bounded(D) || error("complex is not known to be horizontally bounded")
  if horizontal_typ(D) == :chain
    return right_bound(D):-1:left_bound(D)
  elseif horizontal_typ(D) == :cochain
    return left_bound(D):right_bound(D)
  end
  error("typ not recognized")
end

function vertical_range(D::DoubleComplexOfMorphisms) 
  is_vertically_bounded(D) || error("complex is not known to be vertically bounded")
  if vertical_typ(D) == :chain
    return upper_bound(D):-1:lower_bound(D)
  elseif vertical_typ(D) == :cochain
    return lower_bound(D):upper_bound(D)
  end
  error("typ not recognized")
end

function right_bound(D::DoubleComplexOfMorphisms) 
  isdefined(D, :right_bound) || error("complex has no known horizontal upper bound")
  return D.right_bound
end

function left_bound(D::DoubleComplexOfMorphisms) 
  isdefined(D, :left_bound) || error("complex has no known horizontal lower bound")
  return D.left_bound
end

function upper_bound(D::DoubleComplexOfMorphisms) 
  isdefined(D, :upper_bound) || error("complex has no known vertical upper bound")
  return D.upper_bound
end

function lower_bound(D::DoubleComplexOfMorphisms) 
  isdefined(D, :lower_bound) || error("complex has no known vertical lower bound")
  return D.lower_bound
end

extends_right(D::DoubleComplexOfMorphisms) = D.extends_right
extends_left(D::DoubleComplexOfMorphisms) = D.extends_left
extends_up(D::DoubleComplexOfMorphisms) = D.extends_up
extends_down(D::DoubleComplexOfMorphisms) = D.extends_down

### User facing functionality
function getindex(D::DoubleComplexOfMorphisms, t::Tuple)
  return getindex(D, t...)
end

function getindex(D::DoubleComplexOfMorphisms, i::Int, j::Int)
  # legitimize request
  !extends_right(D) && ((has_right_bound(D) &&  i<=right_bound(D)) || error("index out of range"))
  !extends_left(D) && ((has_left_bound(D) && i>=left_bound(D)) || error("index out of range"))
  !extends_up(D) && ((has_upper_bound(D) && j<=upper_bound(D)) || error("index out of range"))
  !extends_down(D) && ((has_lower_bound(D) && j>=lower_bound(D)) || error("index out of range"))

  # load from cache if applicable
  haskey(D.chains, (i, j)) && return D.chains[(i, j)]

  # produce entry otherwise
  new_chain = D.chain_factory(D, i, j)
  D.chains[(i, j)] = new_chain

  # adjust bounds
  has_right_bound(D) && i>right_bound(D) && (D.right_bound = i)
  has_left_bound(D) && i<left_bound(D) && (D.left_bound = i)
  has_upper_bound(D) && j>upper_bound(D) && (D.upper_bound = j)
  has_lower_bound(D) && j<lower_bound(D) && (D.lower_bound = j)
  return new_chain
end

function vertical_map(D::DoubleComplexOfMorphisms, i::Int, j::Int)
  !extends_right(D) && ((has_right_bound(D) && i<=right_bound(D)) || error("index out of range"))
  !extends_left(D) && ((has_left_bound(D) && i>=left_bound(D)) || error("index out of range"))
  if horizontal_typ(D) == :chain
    !extends_up(D) && ((has_upper_bound(D) && j<=upper_bound(D)) || error("index out of range"))
    !extends_down(D) && ((has_lower_bound(D) && j>lower_bound(D)) || error("index out of range"))
  else
    !extends_up(D) && ((has_upper_bound(D) && j<upper_bound(D)) || error("index out of range"))
    !extends_down(D) && ((has_lower_bound(D) && j>=lower_bound(D)) || error("index out of range"))
  end
  haskey(vertical_map_cache(D), (i, j)) && return vertical_map_cache(D)[(i, j)]
  new_map = D.vertical_map_factory(D, i, j)
  vertical_map_cache(D)[(i, j)] = new_map
  return new_map
end

function horizontal_map(D::DoubleComplexOfMorphisms, i::Int, j::Int)
  !extends_up(D) && ((has_upper_bound(D) && j<=upper_bound(D)) || error("index out of range"))
  !extends_down(D) && ((has_lower_bound(D) && j>=lower_bound(D)) || error("index out of range"))
  if horizontal_typ(D) == :chain
    !extends_right(D) && ((has_right_bound(D) && i<=right_bound(D)) || error("index out of range"))
    !extends_left(D) && ((has_left_bound(D) && i>left_bound(D)) || error("index out of range"))
  else
    !extends_right(D) && ((has_right_bound(D) && i<right_bound(D)) || error("index out of range"))
    !extends_left(D) && ((has_left_bound(D) && i>=left_bound(D)) || error("index out of range"))
  end
  haskey(horizontal_map_cache(D), (i, j)) && return horizontal_map_cache(D)[(i, j)]
  new_map = D.horizontal_map_factory(D, i, j)
  horizontal_map_cache(D)[(i, j)] = new_map
  return new_map
end

