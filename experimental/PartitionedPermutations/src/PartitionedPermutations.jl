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
    Perm,
    symmetric_group,
    join,
    compose,
    lower_points,
    upper_points,
    SetPartition,
    subsets,
    Generic,
    size,
    set_partition,
    partition,
    permutation


export PartitionedPermutation
export partitioned_permutation

export partition
export permutation
export length
export adjusted_length
export enumerate_partitioned_permutations
export factorization_partitioned_permutation


include("PartitionedPermutation.jl")
include("PartitionedPermutationProducts.jl")
include("EnumeratePartitionedPermutations.jl")

end

using .PartitionedPermutations

export PartitionedPermutation
export partitioned_permutation

export partition
export permutation
export length
export adjusted_length
export enumerate_partitioned_permutations
export factorization_partitioned_permutation
