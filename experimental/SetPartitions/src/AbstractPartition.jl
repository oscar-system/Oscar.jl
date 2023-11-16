"""
    AbstractPartition

Abstract type for classical-, colored- and spatial partitions.
"""
abstract type AbstractPartition end

function composition(p::T, q::T) where {T <: AbstractPartition}
    composition_loops(p, q)[1]
end

function ⋅(p::T, q::T) where {T <: AbstractPartition}
    composition(p, q)
end

function *(p::T, q::T) where {T <: AbstractPartition}
    composition(p, q)
end

function ⊗(p::T, q::T) where {T <: AbstractPartition}
    tensor_product(p, q)
end

function adjoint(p::AbstractPartition)
    involution(p)
end
