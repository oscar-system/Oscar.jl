### Abstract type and interface for double complexes
abstract type AbsDoubleComplexOfMorphisms{ChainType, MapType} end

### asking for the entries of the complex
getindex(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int) = underlying_double_complex(dc)[i, j]

### asking for horizontal and vertical maps
horizontal_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int) = horizontal_map(underlying_double_complex(dc), i, j)
vertical_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int) = vertical_map(underlying_double_complex(dc), i, j)
vertical_map(dc::AbsDoubleComplexOfMorphisms, t::Tuple) = vertical_map(dc, t...)
horizontal_map(dc::AbsDoubleComplexOfMorphisms, t::Tuple) = horizontal_map(dc, t...)


### asking for the typ of the row and the column complexes
horizontal_typ(dc::AbsDoubleComplexOfMorphisms) = horizontal_typ(underlying_double_complex(dc))
vertical_typ(dc::AbsDoubleComplexOfMorphisms) = vertical_typ(underlying_double_complex(dc))

### asking for known bounds of the row and column complexes
# These are also used to filter legitimate requests to entries and maps, 
# i.e. one can not ask for an entry or a map of a double complex outside 
# of these bounds. If no bound is given in a specific direction, then every 
# request is considered legitimate there. 
#
# Note that at the moment only rectangular shapes are available to describe bounds.
is_horizontally_bounded(dc::AbsDoubleComplexOfMorphisms) = has_horizontal_upper_bound(dc) && has_horizontal_lower_bound(dc)
is_vertically_bounded(dc::AbsDoubleComplexOfMorphisms) = has_vertical_upper_bound(dc) && has_vertical_lower_bound(dc)
is_bounded(dc::AbsDoubleComplexOfMorphisms) = is_horizontally_bounded(dc) && is_vertically_bounded(dc)

has_horizontal_upper_bound(D::AbsDoubleComplexOfMorphisms) = has_horizontal_upper_bound(underlying_double_complex(D))
has_horizontal_lower_bound(D::AbsDoubleComplexOfMorphisms) = has_horizontal_lower_bound(underlying_double_complex(D))
has_vertical_upper_bound(D::AbsDoubleComplexOfMorphisms)   = has_vertical_upper_bound(underlying_double_complex(D))
has_vertical_lower_bound(D::AbsDoubleComplexOfMorphisms)   = has_vertical_lower_bound(underlying_double_complex(D))

is_complete(dc::AbsDoubleComplexOfMorphisms) = is_complete(underlying_double_complex(dc))

# The concrete architecture of double complexes is lazy by default. 
# Hence the constructor needs to be provided with the means to produce 
# the entries of a double complex on request. This is achieved by passing
# certain "factories" to the constructor which carry out this production 
# of the entries and the maps on request. In what follows we specify 
# abstract types for these factories and their interface.

# An abstract type to produce the chains in the (i,j)-th entry of a double
# complex. The interface is formulated for this abstract type, but the 
# user needs to implement a concrete type for any concrete implementation 
# of a double complex.
abstract type ChainFactory{ChainType} end

# Produce the t = (i, j)-th entry which will then be cached.
function (fac::ChainFactory)(dc::AbsDoubleComplexOfMorphisms, t::Tuple)
  return fac(dc, t...)
end

# A dummy placeholder which must be overwritten.
# The first argument will always be the actual double complex itself, 
# so that the body of the function has access to all data already generated 
# and the other functionality available to this double complex. 
function (fac::ChainFactory{ChainType})(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)::ChainType where {ChainType}
  error("production of the ($i, $j)-th chain not implemented")
end

# An abstract type to produce the chain maps going out of the 
# (i, j)-th entry either in the vertical or the horizontal direction.
# The user needs to implement a concrete instance which then knows 
# in particular whether it's supposed to produce the vertical 
# or the horizontal maps (which are then to be cached).
abstract type ChainMorphismFactory{MorphismType} end

function (fac::ChainMorphismFactory)(dc::AbsDoubleComplexOfMorphisms, t1::Tuple)
  return fac(dc, t1...)
end

# A dummy placeholder which must be overwritten; see below.
function (fac::ChainMorphismFactory{MorphismType})(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)::MorphismType where {MorphismType}
  error("could not construct morphism from ($i, $j)")
end

# A minimal concrete type realizing a double complex.
#
# The design is lazy by default. All entries are produced on 
# request and then cached in dictionaries. For the production 
# the user has to provide "factories" in the above sense. 
mutable struct DoubleComplexOfMorphisms{ChainType, MorphismType<:Map} <: AbsDoubleComplexOfMorphisms{ChainType, MorphismType}
  chains::Dict{Tuple{Int, Int}, <:ChainType}
  horizontal_maps::Dict{Tuple{Int, Int}, <:MorphismType}
  vertical_maps::Dict{Tuple{Int, Int}, <:MorphismType}

  # Possible functions to produce new chains and maps in case of an incomplete complex
  chain_factory::ChainFactory{<:ChainType}
  horizontal_map_factory::ChainMorphismFactory{<:MorphismType}
  vertical_map_factory::ChainMorphismFactory{<:MorphismType}

  # Information about the nature of the complex
  is_complete::Bool
  horizontal_typ::Symbol
  vertical_typ::Symbol

  # Information about boundedness and completeness of the complex
  horizontal_upper_bound::Int
  horizontal_lower_bound::Int
  vertical_upper_bound::Int
  vertical_lower_bound::Int

  function DoubleComplexOfMorphisms(
      chain_factory::ChainFactory{ChainType}, 
      horizontal_map_factory::ChainMorphismFactory{MorphismType}, 
      vertical_map_factory::ChainMorphismFactory{MorphismType};
      horizontal_typ::Symbol=:chain,
      vertical_typ::Symbol=:chain,
      horizontal_upper_bound::Union{Int, Nothing} = nothing,
      horizontal_lower_bound::Union{Int, Nothing} = nothing,
      vertical_upper_bound::Union{Int, Nothing} = nothing,
      vertical_lower_bound::Union{Int, Nothing} = nothing,
      is_complete::Bool=false
    ) where {ChainType, MorphismType}
    result = new{ChainType, MorphismType}(Dict{Tuple{Int, Int}, ChainType}(),
                                          Dict{Tuple{Int, Int}, MorphismType}(),
                                          Dict{Tuple{Int, Int}, MorphismType}(),
                                          chain_factory, horizontal_map_factory, 
                                          vertical_map_factory, is_complete, 
                                          horizontal_typ, vertical_typ)
    horizontal_upper_bound !== nothing && (result.horizontal_upper_bound = horizontal_upper_bound)
    horizontal_lower_bound !== nothing && (result.horizontal_lower_bound = horizontal_lower_bound)
    vertical_upper_bound !== nothing && (result.vertical_upper_bound = vertical_upper_bound)
    vertical_lower_bound !== nothing && (result.vertical_lower_bound = vertical_lower_bound)
    return result
  end
end

