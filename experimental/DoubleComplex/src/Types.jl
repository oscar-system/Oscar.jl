### Abstract type and interface for double complexes
abstract type AbsDoubleComplexOfMorphisms{ChainType, MapType} end

### asking for the entries of the complex
getindex(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int) = underlying_double_complex(dc)[i, j]

function has_index(D::AbsDoubleComplexOfMorphisms, t::Tuple)
  return has_index(D, t...)
end

@doc raw"""
    has_index(D::AbsDoubleComplexOfMorphisms, i::Int, j::Int)

Return `true` if the `(i, j)`-th entry of `D` is already known, `false` otherwise.

If the result is `false`, then it might nevertheless still be possible to compute 
`D[i, j]`; use `can_compute_index` for such queries.
"""
function has_index(D::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
  return has_index(underlying_double_complex(D), i, j)
end

function can_compute_index(D::AbsDoubleComplexOfMorphisms, t::Tuple)
  return can_compute_index(D, t...)
end

@doc raw"""
    can_compute_index(D::AbsDoubleComplexOfMorphisms, i::Int, j::Int)

Returns `true` if the entry `D[i, j]` is known or `D` knows how to compute it.
"""
function can_compute_index(D::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
  return can_compute_index(underlying_double_complex(D), i, j)
end

### asking for horizontal and vertical maps
@doc raw"""
    horizontal_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)

Return the morphism ``dc[i, j] → dc[i ± 1, j]`` (the sign depending on the `horizontal_direction` of `dc`). 
"""
horizontal_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int) = horizontal_map(underlying_double_complex(dc), i, j)
horizontal_map(dc::AbsDoubleComplexOfMorphisms, t::Tuple) = horizontal_map(dc, t...)

@doc raw"""
    has_horizontal_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)

Checks whether the double complex `dc` has the horizontal morphism `dc[i, j] → dc[i ± 1, j]`, 
the sign depending on the `horizontal_direction` of `dc`.

If this returns `false` this might just mean that the map has not been computed, yet. 
Use `can_compute_horizontal_map` to learn whether or not this is possible.
"""
has_horizontal_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int) = has_horizontal_map(underlying_double_complex(dc), i, j)
has_horizontal_map(dc::AbsDoubleComplexOfMorphisms, t::Tuple) = has_horizontal_map(dc, t...)

@doc raw"""
    can_compute_horizontal_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)

Returns `true` if `dc` can compute the horizontal morphism `dc[i, j] → dc[i ± 1, j]`, 
the sign depending on the `horizontal_direction` of `dc`, and `false` otherwise.
"""
can_compute_horizontal_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int) = can_compute_horizontal_map(underlying_double_complex(dc), i, j)
can_compute_horizontal_map(dc::AbsDoubleComplexOfMorphisms, t::Tuple) = can_compute_horizontal_map(dc, t...)

@doc raw"""
    vertical_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)

Return the morphism ``dc[i, j] → dc[i, j ± 1]`` (the sign depending on the `vertical_direction` of `dc`). 
"""
vertical_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int) = vertical_map(underlying_double_complex(dc), i, j)

vertical_map(dc::AbsDoubleComplexOfMorphisms, t::Tuple) = vertical_map(dc, t...)

@doc raw"""
    has_vertical_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)

Checks whether the double complex `dc` has the vertical morphism `dc[i, j] → dc[i, j ± 1]`, 
the sign depending on the `vertical_direction` of `dc`.

If this returns `false` this might just mean that the map has not been computed, yet. 
Use `can_compute_vertical_map` to learn whether or not this is possible.
"""
has_vertical_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int) = has_vertical_map(underlying_double_complex(dc), i, j)
has_vertical_map(dc::AbsDoubleComplexOfMorphisms, t::Tuple) = has_vertical_map(dc, t...)

@doc raw"""
    can_compute_vertical_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)

Returns `true` if `dc` can compute the vertical morphism `dc[i, j] → dc[i, j ± 1]`, 
the sign depending on the `vertical_direction` of `dc`, and `false` otherwise.
"""
can_compute_vertical_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int) = can_compute_vertical_map(underlying_double_complex(dc), i, j)
can_compute_vertical_map(dc::AbsDoubleComplexOfMorphisms, t::Tuple) = can_compute_vertical_map(dc, t...)


### asking for the direction of the row and the column complexes
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

Returns `true` if a universal upper bound ``i ≤ B`` for non-zero `D[i, j]` 
is known; `false` otherwise.
"""
has_right_bound(D::AbsDoubleComplexOfMorphisms) = has_right_bound(underlying_double_complex(D))

@doc raw"""
    has_left_bound(D::AbsDoubleComplexOfMorphisms)

Returns `true` if a universal upper bound ``B ≤ i`` for non-zero `D[i, j]` 
is known; `false` otherwise.
"""
has_left_bound(D::AbsDoubleComplexOfMorphisms) = has_left_bound(underlying_double_complex(D))

@doc raw"""
    has_upper_bound(D::AbsDoubleComplexOfMorphisms)

Returns `true` if a universal upper bound ``j ≤ B`` for non-zero `D[i, j]` 
is known; `false` otherwise.
"""
has_upper_bound(D::AbsDoubleComplexOfMorphisms) = has_upper_bound(underlying_double_complex(D))

@doc raw"""
    has_lower_bound(D::AbsDoubleComplexOfMorphisms)

Returns `true` if a universal upper bound ``B ≤ j`` for non-zero `D[i, j]` 
is known; `false` otherwise.
"""
has_lower_bound(D::AbsDoubleComplexOfMorphisms)   = has_lower_bound(underlying_double_complex(D))

@doc raw"""
    right_bound(D::AbsDoubleComplexOfMorphisms)

Returns a bound ``B`` such that `D[i, j]` can be assumed to be zero 
for ``i > B``. Whether or not requests for `D[i, j]` beyond that bound are 
legitimate can be checked using `can_compute_index`.
"""
right_bound(D::AbsDoubleComplexOfMorphisms) = right_bound(underlying_double_complex(D))

@doc raw"""
    left_bound(D::AbsDoubleComplexOfMorphisms)

Returns a bound ``B`` such that `D[i, j]` can be assumed to be zero 
for ``i < B``. Whether or not requests for `D[i, j]` beyond that bound are 
legitimate can be checked using `can_compute_index`.
"""
left_bound(D::AbsDoubleComplexOfMorphisms) = left_bound(underlying_double_complex(D))

@doc raw"""
    upper_bound(D::AbsDoubleComplexOfMorphisms)

Returns a bound ``B`` such that `D[i, j]` can be assumed to be zero 
for ``j > B``. Whether or not requests for `D[i, j]` beyond that bound are 
legitimate can be checked using `can_compute_index`.
"""
upper_bound(D::AbsDoubleComplexOfMorphisms) = upper_bound(underlying_double_complex(D))

@doc raw"""
    lower_bound(D::AbsDoubleComplexOfMorphisms)

Returns a bound ``B`` such that `D[i, j]` can be assumed to be zero 
for ``j < B``. Whether or not requests for `D[i, j]` beyond that bound are 
legitimate can be checked using `can_compute_index`.
"""
lower_bound(D::AbsDoubleComplexOfMorphisms) = lower_bound(underlying_double_complex(D))

@doc raw"""
    is_complete(dc::AbsDoubleComplexOfMorphisms)

Returns `true` if for all indices `(i, j)` with `has_index(dc, i, j) = true` and 
`dc[i, j]` non-zero, the vertex `(i, j)` is lying on an "island" of non-zero entries 
in the grid of the double complex, which is bounded by either zero entries or 
entries for indices `(i', j')` where `can_compute_index(dc, i', j') = false`.
At least one index `dc[i, j]` must be known for this to return `true`.

        ⋮   ⋮   ⋮   ⋮   ⋮   ⋮   ⋮   ⋮   ⋮   ⋮   ⋮   ⋮   ⋮   ⋮   ⋮   ⋮   
        ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   
    … → ? → ? → ? → ? → ? → ? → ? → ? → ? → ? → ? → 0 → 0 → ? → - → - → …
        ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   
    … → ? → ? → 0 → 0 → 0 → ? → ? → 0 → 0 → 0 → 0 → * → * → 0 → - → - → …
        ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   
    … → 0 → 0 → * → * → * → 0 → ? → 0 → ? → 0 → 0 → * → * → * → - → - → …
        ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   
    … → 0 → * → * → * → 0 → ? → ? → 0 → 0 → ? → 0 → * → 0 → * → - → - → …
        ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   
    … → 0 → * → * → * → 0 → ? → 0 → 0 → 0 → 0 → * → * → 0 → 0 → - → - → …
        ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   
    … → ? → 0 → 0 → 0 → ? → ? → ? → ? → ? → ? → 0 → 0 → ? → ? → - → - → …
        ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   
        ⋮   ⋮   ⋮   ⋮   ⋮   ⋮   ⋮   ⋮   ⋮   ⋮   ⋮   ⋮   ⋮   ⋮   ⋮   ⋮   

Example of a pattern of a double complex with `is_complete = true`. 
    `0` : zero entry
    `-` : entry can not be computed (`can_compute_index` returns `false`)
    `*` : non-zero entry which has been computed
    `?` : entry can be computed, but that has not yet been done

!!! note If the double complex has several of the above "islands", then `is_complete` might 
return `true` even though one or more of the "islands" have not yet been uncovered. 
Use this carefully if your full double complex might be separated by zero entries!
"""
function is_complete(D::AbsDoubleComplexOfMorphisms)
  if has_attribute(D, :is_complete) && get_attribute(D, :is_complete)::Bool
    return true
  end
  todo = keys(D.chains)
  isempty(todo) && return false
  for (i, j) in todo
    iszero(D[i, j]) && continue
    has_index(D, i+1, j) || !can_compute_index(D, i+1, j) || return false
    has_index(D, i-1, j) || !can_compute_index(D, i-1, j) || return false
    has_index(D, i, j+1) || !can_compute_index(D, i, j+1) || return false
    has_index(D, i, j-1) || !can_compute_index(D, i, j-1) || return false
  end
  set_attribute!(D, :is_complete, true)
  return true
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

# Test whether the (i, j)-th entry can be computed.
can_compute(fac::ChainFactory, dc::AbsDoubleComplexOfMorphisms, t::Tuple) = can_compute(fac::ChainFactory, dc::AbsDoubleComplexOfMorphisms, t...)

function can_compute(fac::ChainFactory, dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
  error("testing whether the ($i, $j)-th entry can be computed using $fac has not been implemented; see the programmer's documentation on double complexes for details")
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

# Test whether the (i, j)-th entry can be computed.
can_compute(fac::ChainMorphismFactory, dc::AbsDoubleComplexOfMorphisms, t::Tuple) = can_compute(fac::ChainFactory, dc::AbsDoubleComplexOfMorphisms, t...)

function can_compute(fac::ChainMorphismFactory, dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
  error("testing whether the ($i, $j)-th entry can be computed using $fac has not been implemented; see the programmer's documentation on double complexes for details")
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
  right_bound::Union{Nothing, Int}
  left_bound::Union{Nothing, Int}
  upper_bound::Union{Nothing, Int}
  lower_bound::Union{Nothing, Int}

  # caching
  is_complete::Union{Bool, Nothing} # We can not use Bool and isdefined, 
                   # because this leads to errors for fields of type Bool 
                   # with the field being randomly defined or not. 

  function DoubleComplexOfMorphisms(
      chain_factory::ChainFactory{ChainType}, 
      horizontal_map_factory::ChainMorphismFactory{MorphismType}, 
      vertical_map_factory::ChainMorphismFactory{MorphismType};
      horizontal_direction::Symbol=:chain,
      vertical_direction::Symbol=:chain,
      right_bound::Union{Int, Nothing} = nothing,
      left_bound::Union{Int, Nothing} = nothing,
      upper_bound::Union{Int, Nothing} = nothing,
      lower_bound::Union{Int, Nothing} = nothing
    ) where {ChainType, MorphismType}
    result = new{ChainType, MorphismType}(Dict{Tuple{Int, Int}, ChainType}(),
                                          Dict{Tuple{Int, Int}, MorphismType}(),
                                          Dict{Tuple{Int, Int}, MorphismType}(),
                                          chain_factory, horizontal_map_factory, 
                                          vertical_map_factory,
                                          horizontal_direction, vertical_direction)
    result.right_bound = right_bound
    result.left_bound = left_bound
    result.upper_bound = upper_bound
    result.lower_bound = lower_bound
    result.is_complete = nothing # Set value to 'unknown'
    return result
  end
end

