export ActionPolyRing,
       ActionPolyRingElem,
       DifferencePolyRing,
       DifferencePolyRingElem,
       difference_polynomial_ring,
       ndiffs,
       nelemvars

#######################################
#
#  Data types 
#
#######################################

abstract type ActionPolyRing{T} <: Ring end #Perhaps a mutual struct later

abstract type ActionPolyRingElem{T} <: RingElem end

### Difference ###
mutable struct DifferencePolyRing{T} <: ActionPolyRing{T}
  upolyring::AbstractAlgebra.Generic.UniversalPolyRing{T}
  ndiffs::Int
  nelemvars::Int

  DifferencePolyRing{T}(R::Ring, nelemvars::Int, ndiffs::Int) where {T} = new{T}(universal_polynomial_ring(R, nelemvars), ndiffs, nelemvars)
  DifferencePolyRing{T}(R::Ring, S::Vector{Symbol}, ndiffs::Int) where {T} = new{T}(R, ndiffs, length(S), S)
end

mutable struct DifferencePolyRingElem{T} <: ActionPolyRingElem{T}
  coeffs::Vector{DifferencePolyRingElem{T}}
  parent::DifferencePolyRing{T}
  initial::Any
  leader::Any

  DifferencePolyRingElem{T}(dpr::DifferencePolyRing{T}) where {T} = new{T}(DifferencePolyRingElem{T}[], dpr)

end

elem_type(::Type{DifferencePolyRing{T}}) where {T} = DifferencePolyRingElem{T}

parent_type(::Type{DifferencePolyRingElem{T}}) where {T} = DifferencePolyRing{T}

base_ring_type(::Type{<:ActionPolyRing{T}}) where {T} = parent_type(T)

is_domain_type(::Type{<:ActionPolyRingElem{T}}) where {T} = is_domain_type(T)

is_exact_type(::Type{<:ActionPolyRingElem{T}}) where {T} = is_exact_type(T)


#######################################
#
#  Construction 
#
#######################################

##### Algebras #####

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

##### Elements #####

### Difference ###

### Differential ###

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

#######################################
#
#  Auxillary functions 
#
#######################################

