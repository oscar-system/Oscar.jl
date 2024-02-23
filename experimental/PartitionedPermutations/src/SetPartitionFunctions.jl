"""
    number_of_blocks(V::SetPartition)

Return the number of blocks of a `SetPartition`.

# Examples
```jldoctest
julia> number_of_blocks(SetPartition([1, 2, 3], [2, 1, 3, 3]))
3
```
"""
function number_of_blocks(V::SetPartition)
    # obtain one vector describing the partition V
    vec = vcat(V.upper_points, V.lower_points)

    # return the maximum number in vec
    if length(vec) != 0
        return maximum(vec)
    else
        return 0
    end
end

"""
    <=(V::SetPartition, W::SetPartition)

Check if the set partition `V` is dominated by the set partition `W`. This is the case if every block of `V`
is contained in exactly one block of `W`.

# Examples
```jldoctest
julia> SetPartition([1, 1, 2], [1, 2, 3]) <= SetPartition([1, 1, 2], [1, 2, 1])
true
```
"""
function <=(V::SetPartition, W::SetPartition)
    @req size(V) == size(W) "arguments must have the same size"

    # obtain vectors describing the partitions V and W
    V_vec = vcat(V.upper_points, V.lower_points)
    W_vec = vcat(W.upper_points, W.lower_points)

    # introduce a dictionary to store a mapping from the blocks of V to the blocks of W
    block_map = Dict{Int, Int}()

    for (index, block) in enumerate(V_vec)
        # if the block of the index in V has already been mapped to a block of W,
        # check if the mapping is consistent, otherwise add the mapping
        if (haskey(block_map, block) && W_vec[index] != block_map[block])
            return false
        else
            block_map[block] = W_vec[index]
        end
    end
    return true
end

"""
    cycle_partition(p::Perm{Int})

Return the set partition whose blocks are the cycles of the permutation `p`. This set partition has no lower points.

# Examples
```jldoctest
julia> cycle_partition(Perm([2, 1, 3]))
SetPartition([1, 1, 2], Int64[])
```
"""
function cycle_partition(p::Perm{Int})
    cycle_list = collect(cycles(p))
    n = parent(p).n
    partition_vector = zeros(Int64, n)

    for (index, cycle) in enumerate(cycle_list)
        for element in cycle
            partition_vector[element] = index
        end
    end

    return SetPartition(partition_vector, Int64[])
end

"""
    join(V::SetPartition, W::SetPartition)

Return the join of `V` and `W`.

# Examples
```jldoctest
julia> join(SetPartition([1, 2], [2, 3]), SetPartition([1, 2], [1, 3]))
SetPartition([1, 1], [1, 2])
```
"""
function join(V::SetPartition, W::SetPartition)
    @req length(V.upper_points) == length(W.upper_points) "V and W must have the same number of upper points"
    @req length(V.lower_points) == length(W.lower_points) "V and W must have the same number of lower points"

    _V = SetPartition(vcat(V.upper_points, V.lower_points), vcat(V.upper_points, V.lower_points))
    _W = SetPartition(vcat(W.upper_points, W.lower_points), vcat(W.upper_points, W.lower_points))

    join_V_W = compose(_V, _W).upper_points

    return SetPartition(join_V_W[1:length(V.upper_points)], join_V_W[(length(V.upper_points)+1):end])
end
