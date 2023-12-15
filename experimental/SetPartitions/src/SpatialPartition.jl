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
    return SpatialPartition(partition, dim)
end

function hash(p::SpatialPartition, h::UInt)
    return hash(p.partition, hash(p.dimension, h))
end

function ==(p::SpatialPartition, q::SpatialPartition)
    return p.partition == q.partition && p.dimension == q.dimension
end

function deepcopy_internal(p::SpatialPartition, stackdict::IdDict)
    if haskey(stackdict, p)
        return stackdict[p]
    end
    q = SpatialPartition(deepcopy_internal(p.partition, stackdict), 
                         deepcopy_internal(p.dimension, stackdict))
    stackdict[p] = q
    return q
end

function upper_points(p::SpatialPartition)
    return upper_points(p.partition)
end

function lower_points(p::SpatialPartition)
    return lower_points(p.partition)
end

"""
    tensor_product(p::SpatialPartition, q::SpatialPartition)

Return the tensor product of `p` and `q`.
"""
function tensor_product(p::SpatialPartition, q::SpatialPartition)

    @req p.dimension == q.dimension "p and q have different dimensions in tensor product"

    return SpatialPartition(tensor_product(p.partition, q.partition), p.dimension)
end

"""
    involution(p::SpatialPartition)

Return the involution of `p`.
"""
function involution(p::SpatialPartition)
    return SpatialPartition(involution(p.partition), p.dimension)
end

function is_composable(p::SpatialPartition, q::SpatialPartition)
    return p.dimension == q.dimension && is_composable(p.partition, q.partition)
end

"""
    composition_loops(p::SpatialPartition, q::SpatialPartition)

Return the composition of `p` and `q` as well as the number of removed loops.
"""
function composition_loops(p::SpatialPartition, q::SpatialPartition)

    @req is_composable(p, q) "p and q have different dimensions in composition"

    comp_loops = composition_loops(p.partition, q.partition)

    return (SpatialPartition(comp_loops[1], p.dimension), comp_loops[2])

end
