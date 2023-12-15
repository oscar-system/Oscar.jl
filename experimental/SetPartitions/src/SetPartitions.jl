module SetPartitions

import Base: 
    ==, 
    *,
    adjoint,
    copy,
    hash,
    size,
    deepcopy
    
import Oscar:
    ⊗,
    involution,
    tensor_product,
    @req

export ColoredPartition
export SetPartition
export SpatialPartition

export composition
export composition_loops
export construct_category
export is_balanced
export is_noncrossing
export is_pair
export lower_colors
export lower_points
export num_lower_points
export num_upper_points
export print_trace
export rotation
export upper_colors
export upper_points
export vertical_reflection


include("AbstractPartition.jl")
include("Util.jl")
include("SetPartition.jl")
include("ColoredPartition.jl")
include("SpatialPartition.jl")
include("PartitionProperties.jl")
include("GenerateCategory.jl")
end

using .SetPartitions

export ColoredPartition
export SetPartition
export SpatialPartition

export composition
export composition_loops
export construct_category
export is_balanced
export is_noncrossing
export is_pair
export lower_colors
export lower_points
export num_lower_points
export num_upper_points
export print_trace
export rotation
export upper_colors
export upper_points
export vertical_reflection

using Test
using Oscar
include("../test/Partition-test.jl")
include("../test/PartitionProperties-test.jl")
include("../test/GenerateCategory-test.jl")