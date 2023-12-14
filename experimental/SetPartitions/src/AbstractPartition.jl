"""
    AbstractPartition

Abstract type for classical-, colored- and spatial partitions.
"""
abstract type AbstractPartition end

function composition(p::T, q::T) where {T <: AbstractPartition}
    return composition_loops(p, q)[1]
end

function *(p::T, q::T) where {T <: AbstractPartition}
    return composition(p, q)
end

function ⊗(p::T, q::T) where {T <: AbstractPartition}
    return tensor_product(p, q)
end

function adjoint(p::AbstractPartition)
    return involution(p)
end

function num_upper_points(p::AbstractPartition)
	return length(upper_points(p))
end

function num_lower_points(p::AbstractPartition)
	return length(lower_points(p))
end

function size(p::AbstractPartition)
    return num_upper_points(p) + num_lower_points(p)
end