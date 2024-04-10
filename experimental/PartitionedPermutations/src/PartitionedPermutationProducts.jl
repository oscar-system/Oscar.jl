"""
    *(pp_1::PartitionedPermutation, pp_2::PartitionedPermutation)

Return the product of two partitioned permutations as described in [CMS07](@cite).

# Examples
```jldoctest
julia> x = partitioned_permutation(perm(symmetric_group(3), [1, 2, 3]), [1, 2, 3])
PartitionedPermutation((), SetPartition([1, 2, 3], Int64[]))

julia> y = partitioned_permutation(perm(symmetric_group(3), [2, 1, 3]), [1, 1, 3])
PartitionedPermutation((1,2), SetPartition([1, 1, 2], Int64[]))

julia> x*y
PartitionedPermutation((1,2), SetPartition([1, 1, 2], Int64[]))
```
"""
function *(pp_1::PartitionedPermutation, pp_2::PartitionedPermutation)
    # obtain the partitions and permutations from pp_1 and pp_2
    V_1 = partition(pp_1)
    V_2 = partition(pp_2)
    p_1 = permutation(pp_1)
    p_2 = permutation(pp_2) 

    # compute the join of V_1 and V_2, the composition of p_1 and p_2
    W = join(V_1, V_2)
    W_vec = upper_points(W)
    s = p_2 * p_1
    product_pp = PartitionedPermutation(s, W_vec)

    # return the product of pp_1 and pp_2
    if adjusted_length(pp_1) + adjusted_length(pp_2) == adjusted_length(product_pp)
        return product_pp
    else
        return PartitionedPermutation(perm(symmetric_group(length(pp_1)), 1:length(pp_1)), 
            upper_points(cycle_partition(perm(symmetric_group(length(pp_1)), 1:length(pp_1)))))
    end
end

"""
    factor(pp::PartitionedPermutation)

Return the factorization of `pp` in form of a set of 2-tuples.

# Examples
```jldoctest
julia> length(factor(partitioned_permutation(perm(symmetric_group(3), [2, 1, 3]), [1, 1, 2])))
2
```
"""
function factor(pp::PartitionedPermutation)
    size = length(upper_points(partition(pp)))

    product_pairs = Set{Tuple{PartitionedPermutation, PartitionedPermutation}}()
    for pp_1 in enumerate_partitioned_permutations(size)
        for pp_2 in enumerate_partitioned_permutations(size)
            pp_1 * pp_2 == pp && push!(product_pairs, (pp_1, pp_2))
        end
    end
    return product_pairs
end
