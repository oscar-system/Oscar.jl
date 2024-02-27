module PartitionedPermutations

import Base: 
    ==, 
    *,
    adjoint,
    deepcopy,
    deepcopy_internal,
    hash,
    size,
    <=,
    length

import Oscar:
    @req,
    Perm,
    join,
    compose,
    SetPartition,
    subsets,
    cycles,
    Generic,
    size


export PartitionedPermutation

export <=
export cycle_partition
export join
export length
export length2
export enumerate_partitioned_perm
export factorization_partitioned_permutation


include("SetPartitionFunctions.jl")
include("PartitionedPermutation.jl")
include("PartitionedPermutationProducts.jl")
include("EnumeratePartitionedPermutations.jl")

end

using .PartitionedPermutations

export PartitionedPermutation

export <=
export cycle_partition
export join
export length
export length2
export enumerate_partitioned_perm
export factorization_partitioned_permutation
