### Abstract type and interface for double complexes
abstract type AbsDoubleComplexOfMorphisms{ChainType, MapType} end

### asking for the entries of the complex
getindex(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int) = underlying_double_complex(dc)[i, j]

### asking for horizontal and vertical maps
@doc raw"""
    horizontal_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)

Return the morphism ``dc[i, j] → dc[i ± 1, j]`` (the sign depending on the `horizontal_direction` of `dc`). 
"""
horizontal_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int) = horizontal_map(underlying_double_complex(dc), i, j)
@doc raw"""
    vertical_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)

Return the morphism ``dc[i, j] → dc[i, j ± 1]`` (the sign depending on the `vertical_direction` of `dc`). 
"""
vertical_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int) = vertical_map(underlying_double_complex(dc), i, j)

vertical_map(dc::AbsDoubleComplexOfMorphisms, t::Tuple) = vertical_map(dc, t...)
horizontal_map(dc::AbsDoubleComplexOfMorphisms, t::Tuple) = horizontal_map(dc, t...)


### asking for the typ of the row and the column complexes
@doc raw"""
    horizontal_direction(dc::AbsDoubleComplexOfMorphisms)

Return a symbol `:chain` or `:cochain` depending on whether the morphisms of the rows 
of `dc` decrease or increase the (co-)homological index.
"""
horizontal_direction(dc::AbsDoubleComplexOfMorphisms) = horizontal_direction(underlying_double_complex(dc))
@doc raw"""
    vertical_direction(dc::AbsDoubleComplexOfMorphisms)

Return a symbol `:chain` or `:cochain` depending on whether the morphisms of the columns 
of `dc` decrease or increase the (co-)homological index.
"""
vertical_direction(dc::AbsDoubleComplexOfMorphisms) = vertical_direction(underlying_double_complex(dc))

### asking for known bounds of the row and column complexes
# These are also used to filter legitimate requests to entries and maps, 
# i.e. one can not ask for an entry or a map of a double complex outside 
# of these bounds. If no bound is given in a specific direction, then every 
# request is considered legitimate there. 
#
# Note that at the moment only rectangular shapes are available to describe bounds.
is_horizontally_bounded(dc::AbsDoubleComplexOfMorphisms) = has_right_bound(dc) && has_left_bound(dc)
is_vertically_bounded(dc::AbsDoubleComplexOfMorphisms) = has_upper_bound(dc) && has_lower_bound(dc)
is_bounded(dc::AbsDoubleComplexOfMorphisms) = is_horizontally_bounded(dc) && is_vertically_bounded(dc)

@doc raw"""
    has_right_bound(D::AbsDoubleComplexOfMorphisms)

Returns `true` if a universal upper bound ``i ≤ B`` for the `D[i, j]` 
is known; `false` otherwise.
"""
has_right_bound(D::AbsDoubleComplexOfMorphisms) = has_right_bound(underlying_double_complex(D))

@doc raw"""
    has_left_bound(D::AbsDoubleComplexOfMorphisms)

Returns `true` if a universal upper bound ``B ≤ i`` for the `D[i, j]` 
is known; `false` otherwise.
"""
has_left_bound(D::AbsDoubleComplexOfMorphisms) = has_left_bound(underlying_double_complex(D))

@doc raw"""
    has_upper_bound(D::AbsDoubleComplexOfMorphisms)

Returns `true` if a universal upper bound ``j ≤ B`` for the `D[i, j]` 
is known; `false` otherwise.
"""
has_upper_bound(D::AbsDoubleComplexOfMorphisms) = has_upper_bound(underlying_double_complex(D))

@doc raw"""
    has_lower_bound(D::AbsDoubleComplexOfMorphisms)

Returns `true` if a universal upper bound ``B ≤ j`` for the `D[i, j]` 
is known; `false` otherwise.
"""
has_lower_bound(D::AbsDoubleComplexOfMorphisms)   = has_lower_bound(underlying_double_complex(D))

@doc raw"""
    right_bound(D::AbsDoubleComplexOfMorphisms)

Returns a bound ``B`` such that `D[i, j]` is either zero, unknown, or undefined
for ``i > B``. Whether or not requests for `D[i, j]` beyond that bound are 
legitimate depends on the return value of `extends_right(D)`.
"""
right_bound(D::AbsDoubleComplexOfMorphisms) = right_bound(underlying_double_complex(D))

@doc raw"""
    left_bound(D::AbsDoubleComplexOfMorphisms)

Returns a bound ``B`` such that `D[i, j]` is either zero, unknown, or undefined
for ``i < B``. Whether or not requests for `D[i, j]` beyond that bound are 
legitimate depends on the return value of `extends_left(D)`.
"""
left_bound(D::AbsDoubleComplexOfMorphisms) = left_bound(underlying_double_complex(D))

@doc raw"""
    upper_bound(D::AbsDoubleComplexOfMorphisms)

Returns a bound ``B`` such that `D[i, j]` is either zero, unknown, or undefined
for ``j > B``. Whether or not requests for `D[i, j]` beyond that bound are 
legitimate depends on the return value of `extends_up(D)`.
"""
upper_bound(D::AbsDoubleComplexOfMorphisms) = upper_bound(underlying_double_complex(D))

@doc raw"""
    lower_bound(D::AbsDoubleComplexOfMorphisms)

Returns a bound ``B`` such that `D[i, j]` is either zero, unknown, or undefined
for ``j < B``. Whether or not requests for `D[i, j]` beyond that bound are 
legitimate depends on the return value of `extends_down(D)`.
"""
lower_bound(D::AbsDoubleComplexOfMorphisms) = lower_bound(underlying_double_complex(D))

@doc raw"""
    extends_right(D::AbsDoubleComplexOfMorphisms)

Returns `true` if `D` knows how to extend itself to the right, i.e. to produce entries 
`D[i, j]` and (co-)boundary maps for `i > horizontal_right_bound(D)`.
"""
extends_right(D::AbsDoubleComplexOfMorphisms) = extends_right(underlying_double_complex(D))

@doc raw"""
    extends_left(D::AbsDoubleComplexOfMorphisms)

Returns `true` if `D` knows how to extend itself to the left, i.e. to produce entries 
`D[i, j]` and (co-)boundary maps for `i < horizontal_left_bound(D)`.
"""
extends_left(D::AbsDoubleComplexOfMorphisms) = extends_left(underlying_double_complex(D))

@doc raw"""
    extends_up(D::AbsDoubleComplexOfMorphisms)

Returns `true` if `D` knows how to extend itself upwards, i.e. to produce entries 
`D[i, j]` and (co-)boundary maps for `j > upper_bound(D)`.
"""
extends_up(D::AbsDoubleComplexOfMorphisms) = extends_up(underlying_double_complex(D))

@doc raw"""
    extends_down(D::AbsDoubleComplexOfMorphisms)

Returns `true` if `D` knows how to extend itself downwards, i.e. to produce entries 
`D[i, j]` and (co-)boundary maps for `j < lower_bound(D)`.
"""
extends_down(D::AbsDoubleComplexOfMorphisms) = extends_down(underlying_double_complex(D))

@doc raw"""
    is_complete(dc::AbsDoubleComplexOfMorphisms)

Returns `true` if the double complex `dc` has bounds in every direction and 
the double complex can not extend to either direction. 

The entries `dc[i, j]` for `(i, j)` within this range should then be considered 
to comprise all non-zero ones.
"""
function is_complete(dc::AbsDoubleComplexOfMorphisms)
  return has_left_bound(dc) && !extends_left(dc) && has_right_bound(dc) && !extends_right(dc) && has_upper_bound(dc) && !extends_up(dc) && has_lower_bound(dc) && !extends_down(dc)
end

@doc raw"""
    is_horizontally_complete(dc::AbsDoubleComplexOfMorphisms)

Returns `true` if the double complex `dc` has bounds in the horizontal directions and 
it can not extend to either of these directions.

The entries `dc[i, j]` for `(i, j)` within this (possibly unbounded) strip should then be considered 
to comprise all non-zero ones.
"""
function is_horizontally_complete(dc::AbsDoubleComplexOfMorphisms)
  return has_left_bound(dc) && !extends_left(dc) && has_right_bound(dc) && !extends_right(dc) 
end

@doc raw"""
    is_vertically_complete(dc::AbsDoubleComplexOfMorphisms)

Returns `true` if the double complex `dc` has bounds in the vertical directions and 
it can not extend to either of these directions.

The entries `dc[i, j]` for `(i, j)` within this (possibly unbounded) strip should then be considered 
to comprise all non-zero ones.
"""
function is_vertically_complete(dc::AbsDoubleComplexOfMorphisms)
  return has_left_bound(dc) && !extends_left(dc) && has_right_bound(dc) && !extends_right(dc) 
end

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
  horizontal_direction::Symbol
  vertical_direction::Symbol

  # Information about boundedness and completeness of the complex
  right_bound::Int
  left_bound::Int
  upper_bound::Int
  lower_bound::Int

  extends_right::Bool
  extends_left::Bool
  extends_up::Bool
  extends_down::Bool

  function DoubleComplexOfMorphisms(
      chain_factory::ChainFactory{ChainType}, 
      horizontal_map_factory::ChainMorphismFactory{MorphismType}, 
      vertical_map_factory::ChainMorphismFactory{MorphismType};
      horizontal_direction::Symbol=:chain,
      vertical_direction::Symbol=:chain,
      right_bound::Union{Int, Nothing} = nothing,
      left_bound::Union{Int, Nothing} = nothing,
      upper_bound::Union{Int, Nothing} = nothing,
      lower_bound::Union{Int, Nothing} = nothing,
      extends_right::Bool=true,
      extends_left::Bool=true,
      extends_up::Bool=true,
      extends_down::Bool=true
    ) where {ChainType, MorphismType}
    result = new{ChainType, MorphismType}(Dict{Tuple{Int, Int}, ChainType}(),
                                          Dict{Tuple{Int, Int}, MorphismType}(),
                                          Dict{Tuple{Int, Int}, MorphismType}(),
                                          chain_factory, horizontal_map_factory, 
                                          vertical_map_factory,
                                          horizontal_direction, vertical_direction)
    right_bound !== nothing && (result.right_bound = right_bound)
    left_bound !== nothing && (result.left_bound = left_bound)
    upper_bound !== nothing && (result.upper_bound = upper_bound)
    lower_bound !== nothing && (result.lower_bound = lower_bound)
    result.extends_right = extends_right
    result.extends_left = extends_left
    result.extends_up = extends_up
    result.extends_down = extends_down
    return result
  end
end

