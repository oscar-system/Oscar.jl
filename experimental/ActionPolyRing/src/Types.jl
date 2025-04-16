export ActionPolyRing,
       DifferencePolyRing,
       DifferentialPolyRing,
       difference_polynomial_ring,
       differential_polynomial_ring

abstract type ActionPolyRing{T} <: Ring end

mutable struct DifferencePolyRing{T} <: ActionPolyRing{T}
  base_ring::Ring
  nevars::Int64
  S::Vector{Symbol}

  DifferencePolyRing{T}(R::Ring, num_of_elem_vars::Int) where {T} = new{T}(R, num_of_elem_vars, map(x -> Symbol("u" * string(x)), 1:num_of_elem_vars))
  
  DifferencePolyRing{T}(R::Ring, S::Vector{Symbol}) where {T} = new{T}(R, length(S), S)
end

mutable struct DifferentialPolyRing{T} <: ActionPolyRing{T}
  base_ring::Ring
  nevars::Int64
  S::Vector{Symbol}

  DifferentialPolyRing{T}(R::Ring, num_of_elem_vars::Int) where {T} = new{T}(R, num_of_elem_vars, map(x -> Symbol("u" * string(x)), 1:num_of_elem_vars))
  
  DifferentialPolyRing{T}(R::Ring, S::Vector{Symbol}) where {T} = new{T}(R, length(S), S)
end

##### Construction #####

### Difference ###
@doc raw"""
    difference_polynomial_ring(R::Ring, nevars::Int) -> DifferencePolyRing

Return the difference polynomial ring over the ring 'R' in 'nevars' elementary variables.
"""
difference_polynomial_ring(R::Ring, nevars::Int) = DifferencePolyRing{elem_type(R)}(R, nevars)

@doc raw"""
    difference_polynomial_ring(R::Ring, S::Vector{Symbol}) -> DifferencePolyRing

Return the difference polynomial ring over the ring 'R' in $length(S)$ elementary variables with names specified in 'S'.
"""
difference_polynomial_ring(R::Ring, S::Vector{Symbol}) = DifferencePolyRing{elem_type(R)}(R, S)

### Differential ###
@doc raw"""
    differential_polynomial_ring(R::Ring, nevars::Int) -> DifferentialPolyRing

Return the differential polynomial ring over the ring 'R' in 'nevars' elementary variables.
"""
differential_polynomial_ring(R::Ring, nevars::Int) = DifferentialPolyRing{elem_type(R)}(R, nevars)

@doc raw"""
    differential_polynomial_ring(R::Ring, S::Vector{Symbol}) -> DifferentialPolyRing

Return the differential polynomial ring over the ring 'R' in $length(S)$ elementary variables with names specified in 'S'.
"""
differential_polynomial_ring(R::Ring, S::Vector{Symbol}) = DifferentialPolyRing{elem_type(R)}(R, S)


