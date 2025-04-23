export ActionPolyRing,
       ActionPolyRingElem,
       DifferencePolyRing,
       DifferencePolyRingElem,
       DifferentialPolyRing,
       DifferentialPolyRingElem,
       difference_polynomial_ring,
       differential_polynomial_ring,
       ndiffs

#######################################
#
#  Data types 
#
#######################################

abstract type ActionPolyRing{T} <: Ring end #Perhaps a mutual struct later

abstract type ActionPolyRingElem{T} <: RingElem end

### Difference ###
mutable struct DifferencePolyRing{T} <: ActionPolyRing{T}
  base_ring::Ring
  ndiffs::Int
  nelemvars::Int
  S::Vector{Symbol}

  DifferencePolyRing{T}(R::Ring, nelemvars::Int, ndiffs::Int) where {T} = new{T}(R, ndiffs, nelemvars, map(x -> Symbol("u" * string(x)), 1:nelemvars))
  
  DifferencePolyRing{T}(R::Ring, S::Vector{Symbol}, ndiffs::Int) where {T} = new{T}(R, ndiffs, length(S), S)
end

mutable struct DifferencePolyRingElem{T} <: ActionPolyRingElem{T}
  coeffs::Any
  initial::Any
  leader::Symbol
  parent::DifferencePolyRing{T}
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

mutable struct DifferentialPolyRingElem{T} <: ActionPolyRingElem{T}
  coeffs::Any
  initial::Any
  leader::Symbol
  parent::DifferentialPolyRing{T}
end

elem_type(::Type{DifferencePolyRing{T}}) where {T} = DifferencePolyRingElem{T}
elem_type(::Type{DifferentialPolyRing{T}}) where {T} = DifferentialPolyRingElem{T}

parent_type(::Type{DifferencePolyRingElem{T}}) where {T} = DifferencePolyRing{T}
parent_type(::Type{DifferentialPolyRingElem{T}}) where {T} = DifferentialPolyRing{T}

base_ring_type(::Type{DifferencePolyRing{T}}) where {T} = parent_type(T)
base_ring_type(::Type{DifferentialPolyRing{T}}) where {T} = parent_type(T)

is_domain_type(::Type{DifferencePolyRingElem{T}}) where {T} = is_domain_type(T)
is_domain_type(::Type{DifferentialPolyRingElem{T}}) where {T} = is_domain_type(T)

is_exact_type(::Type{DifferencePolyRingElem{T}}) where {T} = is_exact_type(T)
is_exact_type(::Type{DifferentialPolyRingElem{T}}) where {T} = is_exact_type(T)

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
difference_polynomial_ring(R::Ring, nelemvars::Int, ndiffs::Int) = DifferencePolyRing{elem_type(R)}(R, nelemvars, ndiffs)

@doc raw"""
    difference_polynomial_ring(R::Ring, S::Vector{Symbol}, ndiffs::Int) -> DifferencePolyRing

Return the difference polynomial ring over the ring 'R' in $length(S)$ elementary variables with names specified in 'S' and 'ndiffs' commuting endomorphisms.
"""
difference_polynomial_ring(R::Ring, S::Vector{Symbol}, ndiffs::Int) = DifferencePolyRing{elem_type(R)}(R, S, ndiffs)

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

#######################################
#
#  Basic field access 
#
#######################################

### Generic ###
@doc raw"""
    base_ring(apr::ActionPolyRing) -> Ring

Return the base ring of the action polynomial ring 'apr'
"""
base_ring(apr::ActionPolyRing) = apr.base_ring::base_ring_type(typeof(apr))

ndiffs(apr::ActionPolyRing) = apr.ndiffs

nelemvars(apr::ActionPolyRing) = apr.nelemvars

symbols(apr::ActionPolyRing) = apr.S

### Difference ###

### Differential ###
