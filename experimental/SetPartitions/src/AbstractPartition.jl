"""
    AbstractPartition

Abstract type for `SetPartition`, `ColoredPartition` and `SpatialPartition`.
"""
abstract type AbstractPartition end

function compose(p::T, q::T) where {T <: AbstractPartition}
    return compose_count_loops(p, q)[1]
end

function *(p::T, q::T) where {T <: AbstractPartition}
    return compose(p, q)
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