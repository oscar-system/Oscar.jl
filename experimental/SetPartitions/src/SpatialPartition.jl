"""
    SpatialPartition

A spatial partition is a set-partition where the points are 
arranged in multiple layers along the third dimension. 
See Section 2.2 in [CW16](@cite).

It is represented by a set-partition `partition` where the 
number of upper and lower points are a multiple of `dimension`.
"""
struct SpatialPartition <: AbstractPartition
    partition::SetPartition
    dimension::Int
end

function spatial_partition(partition::SetPartition, dim::Int)
    SpatialPartition(partition, dim)
end

function hash(p::SpatialPartition, h::UInt)
    hash(p.partition, hash(p.dimension, h))
end

function ==(p::SpatialPartition, q::SpatialPartition)
    p.partition == q.partition && p.dimension == q.dimension
end

function copy(p::SpatialPartition)
    SpatialPartition(copy(p.partition), copy(p.dimension))
end

"""
    tensor_product(p::SpatialPartition, q::SpatialPartition)

Return the tensor product of `p` and `q`.
"""
function tensor_product(p::SpatialPartition, q::SpatialPartition)

    p.dimension != q.dimension ? 
        error("p and q have different dimensions in tensor product") : 

    SpatialPartition(tensor_product(p.partition, q.partition), p.dimension)
end

"""
    involution(p::SpatialPartition)

Return the involution of `p`.
"""
function involution(p::SpatialPartition)
    SpatialPartition(involution(p.partition), p.dimension)
end

"""
    composition_loops(p::SpatialPartition, q::SpatialPartition)

Return the composition of `p` and `q` as well as the number of removed loops.
"""
function composition_loops(p::SpatialPartition, q::SpatialPartition)

    !is_composable(p, q) ? error("p and q have different dimensions in composition") : 

    comp_loops = composition_loops(p.partition, q.partition)

    (SpatialPartition(comp_loops[1], p.dimension), comp_loops[2])

end
