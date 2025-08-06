"""
    _enumerate_all_partitions(n::Int)

Return and efficiently calculate vector of all partitions of length `n`.
This function is a helper function for `enumerate_partitioned_permutations`.
```
"""
function _enumerate_all_partitions(n::Int)

    partitions = [[1]]
    index_to_cut = 0
    count = false
    
    for i in partitions
        blocks = maximum(i)
        length(i) == n && continue
        count = (length(i) == n - 1)
        for ii in 1:(blocks+1)
            push!(partitions, vcat(i, ii))
            if count
                index_to_cut += 1
            end
        end
    end
    return last(partitions, index_to_cut == 0 ? 1 : index_to_cut)
end

"""
    enumerate_partitioned_permutations(n::Int)

Return and calculate all `PartitionedPermutation` objects of length `n`

# Examples
```jldoctest
julia> length(enumerate_partitioned_permutations(6))
4051
```
"""
function enumerate_partitioned_permutations(n::Int)
    
    partitioned_permutations = []
    # Iterate over all permutations
    for p in symmetric_group(n)
        cycle_part = cycle_partition(p)
        number_of_cycles = number_of_blocks(cycle_part)
        cycle_part_vec = upper_points(cycle_part)


        # Iterate over all partitions dominating p, 
        # those are obtained by merging blocks of the cycle partition of p
        for block_partition in _enumerate_all_partitions(number_of_cycles)
            
            # merge blocks of cycle_part according to block_partition
            part_vec = map(index -> block_partition[index], cycle_part_vec)

            # add the obtained partitioned permutation to the list
            push!(partitioned_permutations, PartitionedPermutation(deepcopy(p), part_vec; check=false))
        end
    end

    return partitioned_permutations
end
