module SetPartitions

import Base: 
    hash, 
    ==, 
    copy

import Oscar:
    involution,
    tensor_product,
    size,
    *,
    ⊗,
    adjoint

export SetPartition
export ColoredPartition
export SpatialPartition
export tensor_product
export composition
export composition_loops
export involution
export rotation
export vertical_reflection
export is_balanced
export is_pair
export is_noncrossing
export construct_category
export print_trace
export size
export *
export ⊗
export ⋅
export adjoint

include("AbstractPartition.jl")
include("Util.jl")
include("SetPartition.jl")
include("ColoredPartition.jl")
include("SpatialPartition.jl")
include("PartitionProperties.jl")
include("GenerateCategory.jl")
end

using .SetPartitions

export SetPartition
export SpatialPartition
export ColoredPartition
export tensor_product
export composition
export composition_loops
export involution
export rotation
export vertical_reflection
export is_balanced
export is_pair
export is_noncrossing
export construct_category
export print_trace
export size
export *
export ⊗
export ⋅
export adjoint
