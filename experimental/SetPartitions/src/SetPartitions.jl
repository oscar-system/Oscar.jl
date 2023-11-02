module SetPartitions

export Partition
export ColoredPartition
export SpatialPartition
export tensor_product
export composition
export involution
export rotation
export vertical_reflection
export is_balanced
export is_pair
export is_noncrossing
export construct_category
export get_trace
export size

include("AbstractPartition.jl")
include("Util.jl")
include("Partition.jl")
include("ColoredPartition.jl")
include("SpatialPartition.jl")
include("PartitionProperties.jl")
include("GenerateCategory.jl")
end

# using .SetPartitions
# 
# export Partition
# export SpatialPartition
# export ColoredPartition
# export tensor_product
# export composition
# export involution
# export rotation
# export vertical_reflection
# export is_balanced
# export is_pair
# export is_noncrossing
# export construct_category
# export get_trace
# export size