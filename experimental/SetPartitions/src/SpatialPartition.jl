import Base.hash
import Base.==
import Base.copy

"""
SpatialPartition

Initialize Spatial Partition object

# Arguments
- `partition`: SetPartition object which is generalized with dimension `dimension`
- `dimension`: dimension of the spatial Partition
"""
struct SpatialPartition <: AbstractPartition
    partition::SetPartition
    dimension::Int
end

function spatial_partition(partition::SetPartition, dim::Int)

    return SpatialPartition(partition, dim)

end

function hash(p::SpatialPartition, h::UInt)

    hash(p.partition, hash(p.dimension, h))
    
end

function ==(p::SpatialPartition, q::SpatialPartition)

    p.partition == q.partition && p.dimension == q.dimension

end

function copy(p::SpatialPartition)
    return SpatialPartition(copy(p.partition), copy(p.dimension))
end

"""
tensor_product(p::SpatialPartition, q::SpatialPartition)

This function applies on p tensor product with q.

# Arguments
- `p`: Input Spatial partition
- `q`: Second input Spatial partition

# Returns
- `p` tensor product `q`
"""
function tensor_product(p::SpatialPartition, q::SpatialPartition)

    p.dimension != q.dimension ? error("p and q have different dimensions in tensor product") : 

    SpatialPartition(tensor_product(p.partition, q.partition), p.dimension)
end

"""
involution(p::SpatialPartition)

This function applies an involution on `p`.

# Arguments
- `p`: Input spatial partition

# Returns
- involution of `p`
"""
function involution(p::SpatialPartition)

    SpatialPartition(involution(p.partition), p.dimension)

end

"""
composition_loops(p::SpatialPartition, q::SpatialPartition)

This function applies composition between p and q.

# Arguments
- `p`: Input partition
- `q`: Second input partition

# Returns
- [`p` composition `q`, number of loops]
"""
function composition_loops(p::SpatialPartition, q::SpatialPartition)

    !is_composable(p, q) ? error("p and q have different dimensions in composition") : 

    comp_loops = composition_loops(p.partition, q.partition)

    (SpatialPartition(comp_loops[1], p.dimension), comp_loops[2])

end
