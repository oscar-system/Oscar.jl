"""
    AbstractPartition

Abstract type for `SetPartition`, `ColoredPartition` and `SpatialPartition`.
"""
abstract type AbstractPartition end

"""
    compose(p::T, q::T) where {T <: AbstractPartition}

Return the composition of `p` and `q` if `p` and `q` are composable.

See `is_composable` and `compose_count_loops`.
"""
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

"""
    number_of_upper_points(p::AbstractPartition)

Return the number of upper points of `p`.
"""
function number_of_upper_points(p::AbstractPartition)
	return length(upper_points(p))
end

"""
    number_of_lower_points(p::AbstractPartition)

Return the number of lower points of `p`.
"""
function number_of_lower_points(p::AbstractPartition)
	return length(lower_points(p))
end

"""
    size(p::AbstractPartition)

Return the total number of points of `p`.
"""
function size(p::AbstractPartition)
    return number_of_upper_points(p) + number_of_lower_points(p)
end
