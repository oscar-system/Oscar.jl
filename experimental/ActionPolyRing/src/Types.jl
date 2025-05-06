export ActionPolyRing,
       ActionPolyRingElem,
       DifferencePolyRing,
       DifferenceVariable,
       DifferencePolyRingElem,
       DifferentialPolyRing,
       DifferentialVariable,
       DifferentialPolyRingElem,
       difference_polynomial_ring,
       difference_variable,
       differential_polynomial_ring,
       differential_variable,
       JetVariable,
       ndiffs,
       nelemvars

#######################################
#
#  Data types 
#
#######################################

abstract type ActionPolyRing{T} <: Ring end #Perhaps a mutual struct later

abstract type ActionPolyRingElem{T} <: RingElem end

abstract type JetVariable{T} end

### Difference ###
mutable struct DifferencePolyRing{T} <: ActionPolyRing{T}
  base_ring::Ring
  ndiffs::Int
  nelemvars::Int
  S::Vector{Symbol}

  DifferencePolyRing{T}(R::Ring, nelemvars::Int, ndiffs::Int) where {T} = new{T}(R, ndiffs, nelemvars, map(x -> Symbol("u" * string(x)), 1:nelemvars))
  DifferencePolyRing{T}(R::Ring, S::Vector{Symbol}, ndiffs::Int) where {T} = new{T}(R, ndiffs, length(S), S)
end

mutable struct DifferenceVariable{T} <: JetVariable{T}
  polyring::DifferencePolyRing{T}
  index::Vector{Int}
  S::Symbol

  function DifferenceVariable{T}(dpr::DifferencePolyRing{T}, elemvar::Symbol, index::Vector{Int}) where {T}
    @req elemvar in symbols(dpr) "Invalid elementary variable"
    @req length(index) == ndiffs(dpr) "Invalid length of index vector"
    @req all(i -> i >= 0, index) "All indices must be nonnegative"
    return new{T}(dpr, index, elemvar)
  end
end

mutable struct DifferencePolyRingElem{T} <: ActionPolyRingElem{T}
  coeffs::Vector{DifferencePolyRingElem{T}}
  #initial::Any
  leader::DifferenceVariable{T}
  parent::DifferencePolyRing{T}

  function DifferencePolyRingElem{T}(ld::DifferenceVariable, cf::Vector{DifferencePolyRingElem{T}}) where {T}
    return _helper(5) 
  end
end


### Differential ###
mutable struct DifferentialPolyRing{T} <: ActionPolyRing{T}
  base_ring::Ring
  ndiffs::Int
  nelemvars::Int
  S::Vector{Symbol}

  DifferentialPolyRing{T}(R::Ring, nelemvars::Int, ndiffs::Int) where {T} = new{T}(R, ndiffs, nelemvars, map(x -> Symbol("u" * string(x)), 1:num_of_elem_vars))
  DifferentialPolyRing{T}(R::Ring, S::Vector{Symbol}, ndiffs::Int) where {T} = new{T}(R, ndiffs, length(S), S)
end

mutable struct DifferentialVariable{T} <: JetVariable{T}
  polyring::DifferentialPolyRing{T}
  index::Vector{Int}
  S::Symbol

  function DifferentialVariable{T}(dpr::DifferentialPolyRing{T}, elemvar::Symbol, index::Vector{Int}) where {T}
    @req elemvar in symbols(dpr) "Invalid elementary variable"
    @req length(index) == ndiffs(dpr) "Invalid length of index vector"
    @req all(i -> i >= 0, index) "All indices must be nonnegative"
    return new{T}(dpr, index, elemvar)
  end
end

mutable struct DifferentialPolyRingElem{T} <: ActionPolyRingElem{T}
  coeffs::Vector{DifferentialPolyRingElem{T}}
  #initial::Any
  leader::DifferentialVariable{T}
  parent::DifferentialPolyRing{T}
end

elem_type(::Type{DifferencePolyRing{T}}) where {T} = DifferencePolyRingElem{T}
elem_type(::Type{DifferentialPolyRing{T}}) where {T} = DifferentialPolyRingElem{T}

parent_type(::Type{DifferencePolyRingElem{T}}) where {T} = DifferencePolyRing{T}
parent_type(::Type{DifferentialPolyRingElem{T}}) where {T} = DifferentialPolyRing{T}

base_ring_type(::Type{<:ActionPolyRing{T}}) where {T} = parent_type(T)

is_domain_type(::Type{<:ActionPolyRingElem{T}}) where {T} = is_domain_type(T)

is_exact_type(::Type{<:ActionPolyRingElem{T}}) where {T} = is_exact_type(T)


#######################################
#
#  Construction 
#
#######################################

### Difference ###
@doc raw"""
    difference_polynomial_ring(R::Ring, nelemvars::Int, ndiffs::Int) -> DifferencePolyRing

Return the difference polynomial ring over the ring 'R' in 'nelemvars' elementary variables and 'ndiffs' commuting endomorphisms.
"""
difference_polynomial_ring(R::Ring, nelemvars::Int, ndiffs::Int) = DifferencePolyRing{elem_type(typeof(R))}(R, nelemvars, ndiffs)

@doc raw"""
    difference_polynomial_ring(R::Ring, S::Vector{Symbol}, ndiffs::Int) -> DifferencePolyRing

Return the difference polynomial ring over the ring 'R' in $length(S)$ elementary variables with names specified in 'S' and 'ndiffs' commuting endomorphisms.
"""
difference_polynomial_ring(R::Ring, S::Vector{Symbol}, ndiffs::Int) = DifferencePolyRing{elem_type(typeof(R))}(R, S, ndiffs)

@doc raw"""
    difference_variable(dpr::DifferencePolyRing, elemvar::Symbol, index::Vector{Int}) -> DifferenceVariable

Create the jet variable of the difference polynomial ring $dpr$ with elementary variable specified by $elemvar$ and index $index$.
"""
difference_variable(dpr::DifferencePolyRing, elemvar::Symbol, index::Vector{Int}) = DifferenceVariable{elem_type(base_ring_type(dpr))}(dpr, elemvar, index)

### Differential ###
@doc raw"""
    differential_polynomial_ring(R::Ring, nelemvars::Int, ndiffs::Int) -> DifferentialPolyRing

Return the differential polynomial ring over the ring 'R' in 'nelemvars' elementary variables and 'ndiffs' commuting derivations.
"""
differential_polynomial_ring(R::Ring, nelemvars::Int, ndiffs::Int) = DifferentialPolyRing{elem_type(R)}(R, nelemvars, ndiffs)

@doc raw"""
    differential_polynomial_ring(R::Ring, S::Vector{Symbol}, ndiffs::Int) -> DifferentialPolyRing

Return the differential polynomial ring over the ring 'R' in $length(S)$ elementary variables with names specified in 'S' and 'ndiffs' commuting derivations.
"""
differential_polynomial_ring(R::Ring, S::Vector{Symbol}, ndiffs::Int) = DifferentialPolyRing{elem_type(R)}(R, S, ndiffs)

@doc raw"""
    differential_variable(dpr::DifferentialPolyRing, elemvar::Symbol, index::Vector{Int}) -> DifferentialVariable

Create the jet variable of the differential polynomial ring $dpr$ with elementary variable specified by $elemvar$ and index $index$.
"""
differential_variable(dpr::DifferentialPolyRing, elemvar::Symbol, index::Vector{Int}) = DifferentialVariable{elem_type(base_ring_type(dpr))}(dpr, elemvar, index)

#######################################
#
#  Basic field access 
#
#######################################

### Difference ###
base_ring(dpr::DifferencePolyRing) = dpr.base_ring::base_ring_type(typeof(dpr))

ndiffs(dpr::DifferencePolyRing) = dpr.ndiffs

nelemvars(dpr::DifferencePolyRing) = dpr.nelemvars

symbols(dpr::DifferencePolyRing) = dpr.S

# Variable
base_ring(dv::DifferenceVariable) = dv.polyring

index(dv::DifferenceVariable) = dv.index

symbols(dv::JetVariable) = dv.S

### Differential ###
base_ring(dpr::DifferentialPolyRing) = dpr.base_ring::base_ring_type(typeof(dpr))

ndiffs(dpr::DifferentialPolyRing) = dpr.ndiffs

nelemvars(dpr::DifferentialPolyRing) = dpr.nelemvars

symbols(dpr::DifferentialPolyRing) = dpr.S

# Variable
base_ring(dv::DifferentialVariable) = dv.polyring

index(dv::DifferentialVariable) = dv.index

symbols(dv::DifferentialVariable) = dv.S

#######################################
#
#  Auxillary functions 
#
#######################################

_helper(n::Int) = "It worked"





