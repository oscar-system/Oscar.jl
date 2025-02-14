module SetPartitions

import Base: 
    +,  
    -,
    *,  
    ==, 
    adjoint,
    deepcopy,
    deepcopy_internal,
    hash,
    join,
    size
    
import Oscar:
    PermGroupElem,
    Ring,
    RingElem,
    RingElement,
    ⊗,
    @req,
    base_ring,
    base_ring_type,
    coefficients,
    compose,
    cycles,
    degree,
    elem_type,
    involution,
    iszero,
    parent,
    parent_type,
    tensor_product

export ColoredPartition
export SetPartition
export SpatialPartition
export LinearPartition

export colored_partition
export compose_count_loops
export construct_category
export cycle_partition
export is_balanced
export is_composable
export is_dominated_by
export is_non_crossing
export is_pair
export join
export levels
export linear_partition
export lower_colors
export lower_points
export number_of_blocks
export number_of_lower_points
export number_of_upper_points
export print_trace
export reflect_vertical
export rotate_top_right
export rotate_top_left
export rotate_bottom_left
export rotate_bottom_right
export set_partition
export spatial_partition
export upper_colors
export upper_points



include("AbstractPartition.jl")
include("Util.jl")
include("SetPartition.jl")
include("ColoredPartition.jl")
include("SpatialPartition.jl")
include("PartitionProperties.jl")
include("GenerateCategory.jl")
include("LinearPartition.jl")
end

using .SetPartitions

export ColoredPartition
export SetPartition
export SpatialPartition
export LinearPartition

export colored_partition
export compose_count_loops
export construct_category
export cycle_partition
export is_balanced
export is_composable
export is_dominated_by
export is_non_crossing
export is_pair
export join
export levels
export linear_partition
export lower_colors
export lower_points
export number_of_blocks
export number_of_lower_points
export number_of_upper_points
export print_trace
export reflect_vertical
export rotate_top_right
export rotate_top_left
export rotate_bottom_left
export rotate_bottom_right
export set_partition
export spatial_partition
export upper_colors
export upper_points
