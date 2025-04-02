########################################################################
# Hyper complexes
#
# The most general type of complexes from which everything else can 
# be derived in principle.
#
# See the more documented double complexes below for details.
########################################################################

### Abstract type and interface for hyper complexes
abstract type AbsHyperComplex{ChainType, MapType} end

### asking for the entries of the complex
getindex(HC::AbsHyperComplex, i::Tuple) = underlying_complex(HC)[i]
getindex(HC::AbsHyperComplex, i::Int...) = getindex(HC, i)
has_index(HC::AbsHyperComplex, i::Tuple) = has_index(underlying_complex(HC), i)
has_index(HC::AbsHyperComplex, i::Int...) = has_index(HC, i)
can_compute_index(HC::AbsHyperComplex, i::Tuple) = can_compute_index(underlying_complex(HC), i)
can_compute_index(HC::AbsHyperComplex, i::Int...) = can_compute_index(HC, i)

### accessing the maps
map(HC::AbsHyperComplex, p::Int, i::Tuple) = map(underlying_complex(HC), p, i)
has_map(HC::AbsHyperComplex, p::Int, i::Tuple) = has_map(underlying_complex(HC), p, i)
can_compute_map(HC::AbsHyperComplex, p::Int, i::Tuple) = can_compute_map(underlying_complex(HC), p, i)

### properties
direction(HC::AbsHyperComplex, p::Int) = direction(underlying_complex(HC), p)
dim(HC::AbsHyperComplex) = dim(underlying_complex(HC))
is_complete(HC::AbsHyperComplex) = is_complete(underlying_complex(HC))
has_upper_bound(HC::AbsHyperComplex, p::Int) = has_upper_bound(underlying_complex(HC), p)
has_lower_bound(HC::AbsHyperComplex, p::Int) = has_lower_bound(underlying_complex(HC), p)
upper_bound(HC::AbsHyperComplex, p::Int) = upper_bound(underlying_complex(HC), p)
lower_bound(HC::AbsHyperComplex, p::Int) = lower_bound(underlying_complex(HC), p)

### internal getters
chain_factory(HC::AbsHyperComplex) = chain_factory(underlying_complex(HC))
map_factory(HC::AbsHyperComplex) = map_factory(underlying_complex(HC))

### type gettes 
chain_type(hc::AbsHyperComplex) = chain_type(typeof(hc))
chain_type(::Type{T}) where {CT, T<:AbsHyperComplex{CT}} = CT

morphism_type(hc::AbsHyperComplex) = morphism_type(typeof(hc))
morphism_type(::Type{T}) where {CT, MT, T<:AbsHyperComplex{CT, MT}} = MT

### factories for the chains
abstract type HyperComplexChainFactory{ChainType} end

function (fac::HyperComplexChainFactory)(HC::AbsHyperComplex, i::Tuple)
  error("production of the $i-th entry not implemented for factory of type $(typeof(fac))")
end

function can_compute(fac::HyperComplexChainFactory, HC::AbsHyperComplex, i::Tuple)
  error("testing whether the $i-th entry can be computed using $fac is not implemented; please overwrite this method")
end

### factories for the maps
abstract type HyperComplexMapFactory{MorphismType} end

function (fac::HyperComplexMapFactory)(HC::AbsHyperComplex, p::Int, i::Tuple)
  error("production of the $i-th map in the $p-th direction not implemented for factory of type $(typeof(fac))")
end

function can_compute(fac::HyperComplexMapFactory, HC::AbsHyperComplex, p::Int, i::Tuple)
  error("testing whether the $i-th map in the $p-th direction can be computed using $fac is not implemented; please overwrite this method")
end

### A minimal concrete type for general hypercomplexes
@attributes mutable struct HyperComplex{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType}
  d::Int 
  chains::Dict{Tuple, <:ChainType}
  morphisms::Dict{Tuple, Dict{Int, <:MorphismType}}

  chain_factory::HyperComplexChainFactory{ChainType}
  map_factory::HyperComplexMapFactory{MorphismType}

  directions::Vector{Symbol}

  upper_bounds::Vector{Union{Int, Nothing}}
  lower_bounds::Vector{Union{Int, Nothing}}

  # fields for caching
  is_complete::Union{Bool, Nothing}
  kernel_cache::Dict{<:Tuple{<:Tuple, Int}, <:Map} # Storing the inclusion map for the kernel
  boundary_cache::Dict{<:Tuple{<:Tuple, Int}, <:Map} # Storing the inclusion map of boundaries
  homology_cache::Dict{<:Tuple{<:Tuple, Int}, <:Map} # Storing the maps from the kernels

  function HyperComplex(
      d::Int,
      chain_factory::HyperComplexChainFactory{ChainType},
      map_factory::HyperComplexMapFactory{MorphismType},
      directions::Vector{Symbol};
      upper_bounds::Vector=[nothing for i in 1:d],
      lower_bounds::Vector=[nothing for i in 1:d]
    ) where {ChainType, MorphismType}
    @assert d >= 0 "can not create negative dimensional hypercomplex"
    chains = Dict{Tuple, ChainType}()
    morphisms = Dict{Tuple, Dict{Int, <:MorphismType}}()
    return new{ChainType, MorphismType}(d, chains, morphisms, 
                                        chain_factory, map_factory, directions, 
                                        Vector{Union{Int, Nothing}}(upper_bounds), 
                                        Vector{Union{Int, Nothing}}(lower_bounds),
                                        nothing
                                       )
  end
end

########################################################################
# Simple complexes derived from hypercomplexes
########################################################################

abstract type AbsSimpleComplex{ChainType, MapType} <:AbsHyperComplex{ChainType, MapType} end

getindex(C::AbsSimpleComplex, i::Int) = C[(i,)]
has_index(C::AbsSimpleComplex, i::Int) = has_index(C, (i,))
can_compute_index(C::AbsSimpleComplex, i::Int) = can_compute_index(C, (i,))
map(C::AbsSimpleComplex, i::Int) = map(C, 1, (i,))
map(C::AbsHyperComplex, i::Int) = (isone(dim(C)) ? map(C, 1, (i,)) : error("complex must be one-dimensional"))
has_map(C::AbsSimpleComplex, i::Int) = has_map(C, 1, (i,))
has_map(C::AbsHyperComplex, i::Int) = (isone(dim(C)) ? has_map(C, 1, (i,)) : error("complex must be one-dimensional"))
can_compute_map(C::AbsSimpleComplex, i::Int) = can_compute_map(C, 1, (i,))
can_compute_map(C::AbsHyperComplex, i::Int) = (isone(dim(C)) ? can_compute_map(C, 1, (i,)) : error("complex must be one-dimensional"))

direction(C::AbsSimpleComplex) = direction(C, 1)
direction(C::AbsHyperComplex) = (isone(dim(C)) ? direction(C, 1) : error("complex must be one-dimensional"))
is_chain_complex(C::AbsSimpleComplex) = direction(C) == :chain
is_cochain_complex(C::AbsSimpleComplex) = !is_chain_complex(C)
has_upper_bound(C::AbsSimpleComplex) = has_upper_bound(C, 1)
has_upper_bound(C::AbsHyperComplex) = (isone(dim(C)) ? has_upper_bound(C, 1) : error("complex must be one-dimensional"))
has_lower_bound(C::AbsSimpleComplex) = has_lower_bound(C, 1)
has_lower_bound(C::AbsHyperComplex) = (isone(dim(C)) ? has_lower_bound(C, 1) : error("complex must be one-dimensional"))
upper_bound(C::AbsSimpleComplex) = upper_bound(C, 1)
upper_bound(C::AbsHyperComplex) = (isone(dim(C)) ? upper_bound(C, 1) : error("complex must be one-dimensional"))
lower_bound(C::AbsSimpleComplex) = lower_bound(C, 1)
lower_bound(C::AbsHyperComplex) = (isone(dim(C)) ? lower_bound(C, 1) : error("complex must be one-dimensional"))
Base.range(C::AbsHyperComplex) = (direction(C) == :chain ? (upper_bound(C):-1:lower_bound(C)) : (lower_bound(C):upper_bound(C)))
map_range(C::AbsHyperComplex) = (direction(C) == :chain ? (upper_bound(C):-1:lower_bound(C)+1) : (lower_bound(C):upper_bound(C)-1))

underlying_complex(C::AbsSimpleComplex) = error("underlying_complex not implemented for $C")

@attributes mutable struct SimpleComplexWrapper{ChainType, MapType} <: AbsSimpleComplex{ChainType, MapType}
  hc::AbsHyperComplex{ChainType, MapType}

  function SimpleComplexWrapper(hc::AbsHyperComplex{ChainType, MapType}) where {ChainType, MapType}
    @assert dim(hc) == 1 "hypercomplex must be one-dimensional"
    return new{ChainType, MapType}(hc)
  end
end

underlying_complex(C::SimpleComplexWrapper) = C.hc

########################################################################
# Double complexes
########################################################################

### Abstract type and interface for double complexes
abstract type AbsDoubleComplexOfMorphisms{ChainType, MapType} <: AbsHyperComplex{ChainType, MapType} end

### Extending generic functionality from hypercomplexes for the double complexes
# Mostly this is providing you with the common syntax.
getindex(C::AbsDoubleComplexOfMorphisms, i::Int, j::Int) = C[(i, j)]

horizontal_range(C::AbsDoubleComplexOfMorphisms) = (direction(C, 1) == :chain ? (left_bound(C):-1:left_bound(C)) : (left_bound(C):right_bound(C)))
vertical_range(C::AbsDoubleComplexOfMorphisms) = (direction(C, 2) == :chain ? (upper_bound(C):-1:lower_bound(C)) : (lower_bound(C):upper_bound(C)))

underlying_complex(C::AbsDoubleComplexOfMorphisms) = error("underlying_complex not implemented for $C")


@doc raw"""
    has_index(D::AbsDoubleComplexOfMorphisms, i::Int, j::Int)

Return `true` if the `(i, j)`-th entry of `D` is already known, `false` otherwise.

If the result is `false`, then it might nevertheless still be possible to compute 
`D[i, j]`; use `can_compute_index` for such queries.
"""
function has_index(D::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
  has_index(D, (i, j))
end

@doc raw"""
    can_compute_index(D::AbsDoubleComplexOfMorphisms, i::Int, j::Int)

Return `true` if the entry `D[i, j]` is known or `D` knows how to compute it.
"""
function can_compute_index(D::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
  can_compute_index(D, (i, j))
end

### asking for horizontal and vertical maps
@doc raw"""
    horizontal_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)

Return the morphism ``dc[i, j] → dc[i ± 1, j]`` (the sign depending on the `horizontal_direction` of `dc`). 
"""
horizontal_map(C::AbsDoubleComplexOfMorphisms, i::Int, j::Int) = map(C, 1, (i, j))
horizontal_map(C::AbsDoubleComplexOfMorphisms, t::Tuple) = map(C, 1, t)

@doc raw"""
    has_horizontal_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)

Check whether the double complex `dc` has the horizontal morphism `dc[i, j] → dc[i ± 1, j]`, 
the sign depending on the `horizontal_direction` of `dc`.

If this returns `false` this might just mean that the map has not been computed, yet. 
Use `can_compute_horizontal_map` to learn whether or not this is possible.
"""
has_horizontal_map(C::AbsDoubleComplexOfMorphisms, i::Int, j::Int) = has_map(C, 1, (i, j))
has_horizontal_map(dc::AbsDoubleComplexOfMorphisms, t::Tuple) = has_map(dc, 1, t)

@doc raw"""
    can_compute_horizontal_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)

Return `true` if `dc` can compute the horizontal morphism `dc[i, j] → dc[i ± 1, j]`,
the sign depending on the `horizontal_direction` of `dc`, and `false` otherwise.
"""
can_compute_horizontal_map(C::AbsDoubleComplexOfMorphisms, i::Int, j::Int) = can_compute_map(C, 1, (i, j))
can_compute_horizontal_map(C::AbsDoubleComplexOfMorphisms, t::Tuple) = can_compute_map(C, 1, t)

@doc raw"""
    vertical_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)

Return the morphism ``dc[i, j] → dc[i, j ± 1]`` (the sign depending on the `vertical_direction` of `dc`).
"""
vertical_map(C::AbsDoubleComplexOfMorphisms, i::Int, j::Int) = map(C, 2, (i, j))
vertical_map(C::AbsDoubleComplexOfMorphisms, t::Tuple) = map(C, 2, t)

@doc raw"""
    has_vertical_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)

Check whether the double complex `dc` has the vertical morphism `dc[i, j] → dc[i, j ± 1]`, 
the sign depending on the `vertical_direction` of `dc`.

If this returns `false` this might just mean that the map has not been computed, yet. 
Use `can_compute_vertical_map` to learn whether or not this is possible.
"""
has_vertical_map(C::AbsDoubleComplexOfMorphisms, i::Int, j::Int) = has_map(C, 2, (i, j))
has_vertical_map(C::AbsDoubleComplexOfMorphisms, t::Tuple) = has_map(C, 2, t)

@doc raw"""
    can_compute_vertical_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)

Return `true` if `dc` can compute the vertical morphism `dc[i, j] → dc[i, j ± 1]`, 
the sign depending on the `vertical_direction` of `dc`, and `false` otherwise.
"""
can_compute_vertical_map(C::AbsDoubleComplexOfMorphisms, i::Int, j::Int) = can_compute_map(C, 2, (i, j))
can_compute_vertical_map(C::AbsDoubleComplexOfMorphisms, t::Tuple) = can_compute_map(C, 2, t)


### asking for the direction of the row and the column complexes
@doc raw"""
    horizontal_direction(dc::AbsDoubleComplexOfMorphisms)

Return a symbol `:chain` or `:cochain` depending on whether the morphisms of the rows 
of `dc` decrease or increase the (co-)homological index.
"""
horizontal_direction(dc::AbsDoubleComplexOfMorphisms) = direction(dc, 1)

@doc raw"""
    vertical_direction(dc::AbsDoubleComplexOfMorphisms)

Return a symbol `:chain` or `:cochain` depending on whether the morphisms of the columns 
of `dc` decrease or increase the (co-)homological index.
"""
vertical_direction(dc::AbsDoubleComplexOfMorphisms) = direction(dc, 2)

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

Return `true` if a universal upper bound ``i ≤ B`` for non-zero `D[i, j]`
is known; `false` otherwise.
"""
has_right_bound(C::AbsDoubleComplexOfMorphisms) = has_upper_bound(C, 1)

@doc raw"""
    has_left_bound(D::AbsDoubleComplexOfMorphisms)

Return `true` if a universal upper bound ``B ≤ i`` for non-zero `D[i, j]`
is known; `false` otherwise.
"""
has_left_bound(C::AbsDoubleComplexOfMorphisms) = has_lower_bound(C, 1)

@doc raw"""
    has_upper_bound(D::AbsDoubleComplexOfMorphisms)

Return `true` if a universal upper bound ``j ≤ B`` for non-zero `D[i, j]`
is known; `false` otherwise.
"""
has_upper_bound(C::AbsDoubleComplexOfMorphisms) = has_upper_bound(C, 2)

@doc raw"""
    has_lower_bound(D::AbsDoubleComplexOfMorphisms)

Return `true` if a universal upper bound ``B ≤ j`` for non-zero `D[i, j]`
is known; `false` otherwise.
"""
has_lower_bound(C::AbsDoubleComplexOfMorphisms) = has_lower_bound(C, 2)

@doc raw"""
    right_bound(D::AbsDoubleComplexOfMorphisms)

Return a bound ``B`` such that `D[i, j]` can be assumed to be zero
for ``i > B``. Whether or not requests for `D[i, j]` beyond that bound are 
legitimate can be checked using `can_compute_index`.
"""
right_bound(C::AbsDoubleComplexOfMorphisms) = upper_bound(C, 1)

@doc raw"""
    left_bound(D::AbsDoubleComplexOfMorphisms)

Return a bound ``B`` such that `D[i, j]` can be assumed to be zero 
for ``i < B``. Whether or not requests for `D[i, j]` beyond that bound are 
legitimate can be checked using `can_compute_index`.
"""
left_bound(C::AbsDoubleComplexOfMorphisms) = lower_bound(C, 1)

@doc raw"""
    upper_bound(D::AbsDoubleComplexOfMorphisms)

Return a bound ``B`` such that `D[i, j]` can be assumed to be zero 
for ``j > B``. Whether or not requests for `D[i, j]` beyond that bound are 
legitimate can be checked using `can_compute_index`.
"""
upper_bound(C::AbsDoubleComplexOfMorphisms) = upper_bound(C, 2)

@doc raw"""
    lower_bound(D::AbsDoubleComplexOfMorphisms)

Return a bound ``B`` such that `D[i, j]` can be assumed to be zero
for ``j < B``. Whether or not requests for `D[i, j]` beyond that bound are 
legitimate can be checked using `can_compute_index`.
"""
lower_bound(C::AbsDoubleComplexOfMorphisms) = lower_bound(C, 2)

@doc raw"""
    is_complete(dc::AbsDoubleComplexOfMorphisms)

Return `true` if for all indices `(i, j)` with `has_index(dc, i, j) = true` and
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
  return is_complete(underlying_complex(D))
end


### A concrete type using hypercomplexes in the background
@attributes mutable struct DoubleComplexWrapper{ChainType, MapType} <: AbsDoubleComplexOfMorphisms{ChainType, MapType}
  hc::AbsHyperComplex{ChainType, MapType}

  function DoubleComplexWrapper(hc::AbsHyperComplex{ChainType, MapType}) where {ChainType, MapType}
    @assert dim(hc) == 2 "hypercomplex must be one-dimensional"
    return new{ChainType, MapType}(hc)
  end
end

underlying_complex(dc::DoubleComplexWrapper) = dc.hc

########################################################################
# A standalone concrete type for double complexes 
########################################################################

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
function (fac::ChainFactory{ChainType})(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int) where {ChainType}
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
function (fac::ChainMorphismFactory{MorphismType})(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int) where {MorphismType}
  error("could not construct morphism from ($i, $j)")
end

# Test whether the (i, j)-th entry can be computed.
can_compute(fac::ChainMorphismFactory, dc::AbsDoubleComplexOfMorphisms, t::Tuple) = can_compute(fac::ChainFactory, dc::AbsDoubleComplexOfMorphisms, t...)

function can_compute(fac::ChainMorphismFactory, dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
  error("testing whether the ($i, $j)-th entry can be computed using $fac has not been implemented; see the programmer's documentation on double complexes for details")
end

### The actual concrete type
@attributes mutable struct DoubleComplexOfMorphisms{ChainType, MorphismType} <: AbsDoubleComplexOfMorphisms{ChainType, MorphismType}
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


### Zero dimensional hypercomplexes
# Needed for internals.

struct DummyChainFactory{ChainType} <: HyperComplexChainFactory{ChainType} end
struct DummyMorphismFactory{MorphismType} <: HyperComplexMapFactory{MorphismType} end

@attributes mutable struct ZeroDimensionalComplex{ChainType, MorphismType} <: AbsDoubleComplexOfMorphisms{ChainType, MorphismType}
  M::ChainType

  function ZeroDimensionalComplex(M::ChainType) where {ChainType}
    MorphismType = morphism_type(ChainType)
    return new{ChainType, MorphismType}(M)
  end
end

function getindex(c::ZeroDimensionalComplex, i::Tuple) 
  @assert iszero(length(i)) "index out of bounds"
  return c.M
end

function can_compute_index(c::ZeroDimensionalComplex, i::Tuple)
  @assert iszero(length(i)) "index out of bounds"
  return true
end

function has_index(c::ZeroDimensionalComplex, i::Tuple)
  @assert iszero(length(i)) "index out of bounds"
  return true
end

dim(c::ZeroDimensionalComplex) = 0
is_complete(c::ZeroDimensionalComplex) = true

### Wrappers for shifted complexes
@attributes mutable struct ShiftedHyperComplex{ChainType, MorphismType} <: AbsDoubleComplexOfMorphisms{ChainType, MorphismType}
  complex::AbsHyperComplex{ChainType, MorphismType}
  shift::Vector{Int}

  function ShiftedHyperComplex(c::AbsHyperComplex{ChainType, MorphismType}, shift::Tuple) where {ChainType, MorphismType}
    @assert length(shift) == dim(c)
    return new{ChainType, MorphismType}(c, collect(shift))
  end
end

### Views of hypercomplexes
#
# These can be slices, hyperrectangles, extensions, ...
@attributes mutable struct HyperComplexView{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType}
  d::Int
  original::AbsHyperComplex # the complex from which the view is taken; use `original_complex`
  mapping_matrix::Matrix{Int} # an integer matrix to map an index of the new complex to the 
                              # corresponding index of the old; use `mapping_matrix`
  offset_vector::Vector{Int}  # a vector to be added to the result of the mapping_matrix 
                              # multiplication to arrive at the actual index of the old complex
  lower_bounds::Vector{Int}
  upper_bounds::Vector{Int}

  # v is a vector with entry (k, u) in place i signaling that 
  #   * if u is a UnitRange: the i-th dimension of the 
  #     old complex should become the k-th dimension of 
  #     the new one with index running within u.
  #
  #   * if u is an integer: the i-th dimension of the old 
  #     complex does not appear in the new complex, and 
  #     its value is fixed to be u.
  #
  #   * if i is bigger than the dimension of the old complex:
  #     k is a new dimension of the new complex and the new 
  #     complex is concentrated in some degree u = r:r there.
  function HyperComplexView(
      c::AbsHyperComplex{ChainType, MapType}, v::Vector{T}
    ) where {ChainType, MapType, T<:Tuple{Int, <:Union{Int, UnitRange{Int}}}}
    @assert length(v) >= dim(c) "not enough information provided for a full view"
    @assert all(k->k[2] isa UnitRange && isone(length(k[2])), v[dim(c)+1:end]) "additional dimensions must have a specified range of length one"
    ranges = [(k, u) for (k, u) in v if u isa UnitRange]
    d = length(ranges)
    all_ind = [k for (k, u) in v if u isa UnitRange]
    @assert all(in(all_ind), 1:d) "matching of ranges is not unique"

    mapping_matrix = [0 for i in 1:dim(c), j in 1:d]
    new_ranges = Dict(v)
    offset_vector = [0 for i in 1:dim(c)]
    for i in 1:dim(c)
      (k, u) = v[i]
      if u isa UnitRange
        mapping_matrix[i, k] = 1
      else
        iszero(k) || error("fixed slice indices must be provided without a new dimension")
        offset_vector[i] = u
      end
    end
    lower_bounds = [first(new_ranges[k]) for k in 1:d]
    upper_bounds = [last(new_ranges[k]) for k in 1:d]
    
    return new{ChainType, MapType}(d, c, mapping_matrix, offset_vector, lower_bounds, upper_bounds)
  end
end

@attributes mutable struct TotalComplex{ChainType, MorphismType, ComplexType} <: AbsSimpleComplex{ChainType, MorphismType}
  original::ComplexType
  complex::AbsHyperComplex{ChainType, MorphismType}

  kernels::Dict{Int, ChainType}
  boundaries::Dict{Int, ChainType}
  homologies::Dict{Int, ChainType}

  function TotalComplex(c::ComplexType, internal_complex::AbsHyperComplex{ChainType, MorphismType}) where {ChainType, MorphismType, ComplexType}
    all(k->has_upper_bound(c, k), 1:dim(c)) || all(k->has_lower_bound(c, k), 1:dim(c)) || error("complex must have either all upper or all lower bounds")
    return new{ChainType, MorphismType, ComplexType}(c, internal_complex)
  end
end

underlying_complex(c::TotalComplex) = c.complex

