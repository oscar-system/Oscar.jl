####################################################################
# Functions for Set Partitions
####################################################################
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



############################################################
# Partitioned Permutations
############################################################

"""
    PartitionedPermutation

The type of partitioned permutations. Fieldnames are
- p::Perm{Int} - a permutation
- V::SetPartition - a partition
- check::Bool = true
If the permutation has length `n`, then the partition must have `n` upper points and 0 lower points. 
Further, if `W` is the partition given by the cycles of `p`, then `W` must be dominated by `V` in the 
sense that every block of `W` is contained in one block of `V`. There is one inner constructer of PartitionedPermutation:
- PartitionedPermutation(_p::Perm{Int}, _V::Vector{Int}) constructs the partitioned permutation where the partition is given by the vector _V.
If the optional flag `check` is set to `false`, then the constructor skips the validation of the requirements mentioned above.

# Examples
```jldoctest
julia> PartitionedPermutation(Perm([2, 1, 3]), [1, 1, 2])
PartitionedPermutation((1,2), SetPartition([1, 1, 2], Int64[]))
```
"""
struct PartitionedPermutation
    p::Perm{Int}
    V::SetPartition
    check::Bool

    function PartitionedPermutation(_p::Perm{Int}, _V::Vector{Int}, _check::Bool=true)
        __V = SetPartition(_V, Int[])
        if _check
            @req parent(_p).n == length(_V) "permutation and partition must have the same length"
            @req cycle_partition(_p) <= __V "permutation must be dominated by partition"
        end
        new(_p, __V)
    end
end

function partitioned_permutation(p::Perm{Int}, V::Vector{Int}, check::Bool=true)
    return PartitionedPermutation(p, V, check)
end

function ==(pp_1::PartitionedPermutation, pp_2::PartitionedPermutation)
    return pp_1.p == pp_2.p && pp_1.V == pp_2.V
end

function hash(pp::PartitionedPermutation, h::UInt)
    return hash(pp.p, hash(pp.V, h))
end

function deepcopy_internal(pp::PartitionedPermutation, stackdict::IdDict)
    if haskey(stackdict, pp)
        return stackdict[pp]
    end
    q = PartitionedPermutation(deepcopy_internal(pp.p, stackdict), 
                     deepcopy_internal(pp.V, stackdict))
    stackdict[pp] = q
    return q
end

"""
    length(pp::PartitionedPermutation)

Return the length of a partitioned permutation, i.e. the size of the underlying set.

# Examples
```jldoctest
julia> length(PartitionedPermutation(Perm([2, 1]), [1, 1]))
2
```
"""
function length(pp::PartitionedPermutation)
    return parent(pp.p).n
end

"""
    length2(pp::PartitionedPermutation)

Return the adjusted length of a partitioned permutation as described in [CITE] as `|(V, pi)|`
for a partition `V` and a permutation `pi`.

# Examples
```jldoctest
julia> length2(PartitionedPermutation(Perm([2, 1]), [1, 1]))
1
```
"""
function length2(pp::PartitionedPermutation)
    p = pp.p
    V = pp.V
    return parent(p).n - (2*number_of_blocks(V) - length(cycles(p)))
end

"""
    *(pp_1::PartitionedPermutation, pp_2::PartitionedPermutation)

Return the product of two partitioned permutations as described in [CITE].

# Examples
```jldoctest
julia> x = PartitionedPermutation(Perm([1, 2, 3]), [1, 2, 3])
PartitionedPermutation((), SetPartition([1, 2, 3], Int64[]))

julia> y = PartitionedPermutation(Perm([2, 1, 3]), [1, 1, 3])
PartitionedPermutation((1,2), SetPartition([1, 1, 2], Int64[]))

julia> x*y
PartitionedPermutation((1,2), SetPartition([1, 1, 2], Int64[]))
```
"""
function *(pp_1::PartitionedPermutation, pp_2::PartitionedPermutation)
    # obtain the partitions and permutations from pp_1 and pp_2
    V_1 = pp_1.V 
    V_2 = pp_2.V
    p_1 = pp_1.p
    p_2 = pp_2.p 

    # compute the join of V_1 and V_2, the composition of p_1 and p_2
    W = join(V_1, V_2)
    W_vec = W.upper_points
    s = p_2 * p_1
    product_pp = PartitionedPermutation(s, W_vec)

    # return the product of pp_1 and pp_2
    if length2(pp_1) + length2(pp_2) == length2(product_pp)
        return product_pp
    else
        return PartitionedPermutation(Perm(1:length(pp_1)), cycle_partition(Perm(1:length(pp_1))).upper_points)
    end
end
