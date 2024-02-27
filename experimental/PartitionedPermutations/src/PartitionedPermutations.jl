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
    Perm,
    symmetric_group,
    join,
    compose,
    lower_points,
    upper_points,
    SetPartition,
    subsets,
    cycles,
    Generic,
    size


export PartitionedPermutation

export is_dominated_by
export cycle_partition
export join
export length
export length2
export enumerate_partitioned_permutations
export factorization_partitioned_permutation


include("SetPartitionFunctions.jl")
include("PartitionedPermutation.jl")
include("PartitionedPermutationProducts.jl")
include("EnumeratePartitionedPermutations.jl")

end

using .PartitionedPermutations

export PartitionedPermutation

export is_dominated_by
export cycle_partition
export join
export length
export length2
export enumerate_partitioned_permutations
export factorization_partitioned_permutation
