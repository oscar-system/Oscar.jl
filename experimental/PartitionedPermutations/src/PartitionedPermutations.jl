module PartitionedPermutations

import Base: 
    ==, 
    *,
    adjoint,
    deepcopy,
    deepcopy_internal,
    hash,
    size,
    length

import Oscar:
    @req,
    cycle_partition,
    cycles,
    is_dominated_by,
    number_of_blocks,
    PermGroupElem,
    perm,
    symmetric_group,
    degree,
    join,
    compose,
    lower_points,
    upper_points,
    SetPartition,
    subsets,
    size,
    set_partition,
    partition,
    permutation,
    factor


export PartitionedPermutation
export partitioned_permutation

export partition
export permutation
export length
export adjusted_length
export enumerate_partitioned_permutations


include("PartitionedPermutation.jl")
include("PartitionedPermutationProducts.jl")
include("EnumeratePartitionedPermutations.jl")

end

using .PartitionedPermutations

export PartitionedPermutation
export partitioned_permutation

export adjusted_length
export enumerate_partitioned_permutations
