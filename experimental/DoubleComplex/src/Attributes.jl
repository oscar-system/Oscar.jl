### Basic getters
function is_complete(D::DoubleComplexOfMorphisms)
  return D.is_complete
end

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
has_horizontal_upper_bound(D::DoubleComplexOfMorphisms) = isdefined(D, :horizontal_upper_bound)
has_horizontal_lower_bound(D::DoubleComplexOfMorphisms) = isdefined(D, :horizontal_lower_bound)
has_vertical_upper_bound(D::DoubleComplexOfMorphisms)   = isdefined(D, :vertical_upper_bound)
has_vertical_lower_bound(D::DoubleComplexOfMorphisms)   = isdefined(D, :vertical_lower_bound)

function horizontal_range(D::DoubleComplexOfMorphisms) 
  is_horizontally_bounded(D) || error("complex is not known to be horizontally bounded")
  if horizontal_typ(D) == :chain
    return horizontal_upper_bound(D):-1:horizontal_lower_bound(D)
  elseif horizontal_typ(D) == :cochain
    return horizontal_lower_bound(D):horizontal_upper_bound(D)
  end
  error("typ not recognized")
end

function vertical_range(D::DoubleComplexOfMorphisms) 
  is_vertically_bounded(D) || error("complex is not known to be vertically bounded")
  if vertical_typ(D) == :chain
    return vertical_upper_bound(D):-1:vertical_lower_bound(D)
  elseif vertical_typ(D) == :cochain
    return vertical_lower_bound(D):vertical_upper_bound(D)
  end
  error("typ not recognized")
end

function horizontal_upper_bound(D::DoubleComplexOfMorphisms) 
  isdefined(D, :horizontal_upper_bound) || error("complex has no known horizontal upper bound")
  return D.horizontal_upper_bound
end

function horizontal_lower_bound(D::DoubleComplexOfMorphisms) 
  isdefined(D, :horizontal_lower_bound) || error("complex has no known horizontal lower bound")
  return D.horizontal_lower_bound
end

function vertical_upper_bound(D::DoubleComplexOfMorphisms) 
  isdefined(D, :vertical_upper_bound) || error("complex has no known vertical upper bound")
  return D.vertical_upper_bound
end

function vertical_lower_bound(D::DoubleComplexOfMorphisms) 
  isdefined(D, :vertical_lower_bound) || error("complex has no known vertical lower bound")
  return D.vertical_lower_bound
end

### User facing functionality
function getindex(D::DoubleComplexOfMorphisms, t::Tuple)
  return getindex(D, t...)
end

function getindex(D::DoubleComplexOfMorphisms, i::Int, j::Int)
  (has_horizontal_upper_bound(D) && i<=horizontal_upper_bound(D)) || error("index out of range")
  (has_horizontal_lower_bound(D) && i>=horizontal_lower_bound(D)) || error("index out of range")
  (has_vertical_upper_bound(D) && j<=vertical_upper_bound(D)) || error("index out of range")
  (has_vertical_lower_bound(D) && j>=vertical_lower_bound(D)) || error("index out of range")
  haskey(D.chains, (i, j)) && return D.chains[(i, j)]
  new_chain = D.chain_factory(D, i, j)
  D.chains[(i, j)] = new_chain
  return new_chain
end

function vertical_map(D::DoubleComplexOfMorphisms, i::Int, j::Int)
  (has_horizontal_upper_bound(D) && i<=horizontal_upper_bound(D)) || error("index out of range")
  (has_horizontal_lower_bound(D) && i>=horizontal_lower_bound(D)) || error("index out of range")
  if horizontal_typ(D) == :chain
    (has_vertical_upper_bound(D) && j<=vertical_upper_bound(D)) || error("index out of range")
    (has_vertical_lower_bound(D) && j>vertical_lower_bound(D)) || error("index out of range")
  else
    (has_vertical_upper_bound(D) && j<vertical_upper_bound(D)) || error("index out of range")
    (has_vertical_lower_bound(D) && j>=vertical_lower_bound(D)) || error("index out of range")
  end
  haskey(vertical_map_cache(D), (i, j)) && return vertical_map_cache(D)[(i, j)]
  new_map = D.vertical_map_factory(D, i, j)
  vertical_map_cache(D)[(i, j)] = new_map
  return new_map
end

function horizontal_map(D::DoubleComplexOfMorphisms, i::Int, j::Int)
  (has_vertical_upper_bound(D) && j<=vertical_upper_bound(D)) || error("index out of range")
  (has_vertical_lower_bound(D) && j>=vertical_lower_bound(D)) || error("index out of range")
  if horizontal_typ(D) == :chain
    (has_horizontal_upper_bound(D) && i<=horizontal_upper_bound(D)) || error("index out of range")
    (has_horizontal_lower_bound(D) && i>horizontal_lower_bound(D)) || error("index out of range")
  else
    (has_horizontal_upper_bound(D) && i<horizontal_upper_bound(D)) || error("index out of range")
    (has_horizontal_lower_bound(D) && i>=horizontal_lower_bound(D)) || error("index out of range")
  end
  haskey(horizontal_map_cache(D), (i, j)) && return horizontal_map_cache(D)[(i, j)]
  new_map = D.horizontal_map_factory(D, i, j)
  horizontal_map_cache(D)[(i, j)] = new_map
  return new_map
end

