"""
    AbstractPartition

Abstract type for classical-, colored- and spatial partitions.
"""
abstract type AbstractPartition end

function composition(p::T, q::T) where {T <: AbstractPartition}
    return composition_loops(p, q)[1]
end

function ⋅(p::T, q::T) where {T <: AbstractPartition}
    return composition(p, q)
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
