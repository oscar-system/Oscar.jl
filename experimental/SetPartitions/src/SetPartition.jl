"""
    SetPartition

`SetPartition` represents a partition of a set of upper and lower points 
into disjoint subsets. 
Such set-partitions are often depicted as string diagrams where points in the same subset 
are connected by lines. See also Section 4.1.1 in [Gro20](@cite).
"""
struct SetPartition <: AbstractPartition
    upper_points::Vector{Int}
    lower_points::Vector{Int}

    function SetPartition(upper_points::Vector{Int}, lower_points::Vector{Int})
        (new_upper, new_lower) = _normal_form(upper_points, lower_points)
        return new(new_upper, new_lower)
    end
end

"""
    set_partition(upper_points::Vector, lower_points::Vector)

Construct a `SetPartition` with points `upper_points` and `lower_points`
where two points are in the same subset if and only if they have the same value.

Note that `upper_points` and `lower_points` are stored in normal form.

# Examples
```jldoctest
julia> set_partition([2, 4], [4, 99])
SetPartition([1, 2], [2, 3])
```
"""
function set_partition(upper_points::Vector, lower_points::Vector)
    return SetPartition(convert(Vector{Int}, upper_points), 
                        convert(Vector{Int}, lower_points))
end


function hash(p::SetPartition, h::UInt)
    return hash(upper_points(p), hash(lower_points(p), h))
end

function ==(p::SetPartition, q::SetPartition)
    return upper_points(p) == upper_points(q) && 
           lower_points(p) == lower_points(q) 
    
end

function deepcopy_internal(p::SetPartition, stackdict::IdDict)
    if haskey(stackdict, p)
        return stackdict[p]
    end
    q = SetPartition(deepcopy_internal(p.upper_points, stackdict), 
                     deepcopy_internal(p.lower_points, stackdict))
    stackdict[p] = q
    return q
end

"""
    upper_points(p::SetPartition)

Return the upper points of `p`.

# Examples
```jldoctest
julia> upper_points(set_partition([2, 4], [4, 99]))
2-element Vector{Int64}:
 1
 2
```
"""
function upper_points(p::SetPartition)
    return p.upper_points
end

"""
    lower_points(p::SetPartition)

Return the lower points of `p`.

# Examples
```jldoctest
julia> lower_points(set_partition([2, 4], [4, 99]))
2-element Vector{Int64}:
 2
 3
```
"""
function lower_points(p::SetPartition)
    return p.lower_points
end

"""
    tensor_product(p::SetPartition, q::SetPartition)

Return the tensor product of `p` and `q`.

The tensor product of two partitions is given by their horizontal concatenation.
See also Section 4.1.1 in [Gro20](@cite).

# Examples
```jldoctest
julia> tensor_product(set_partition([1, 2], [2, 1]), set_partition([1, 1], [1]))
SetPartition([1, 2, 3, 3], [2, 1, 3])
```
"""
function tensor_product(p::SetPartition, q::SetPartition)
    q_new = _new_point_values(upper_points(p), lower_points(p), 
                              upper_points(q), lower_points(q))

    return set_partition(vcat(upper_points(p), q_new[1]), 
                         vcat(lower_points(p), q_new[2]))
end

"""
    involution(p::SetPartition)

Return the involution of `p`.

The involution of a partition is obtained by swapping the upper and lower points.
See also Section 4.1.1 in [Gro20](@cite).

# Examples
```jldoctest
julia> involution(set_partition([1, 2, 3], [2, 1]))
SetPartition([1, 2], [2, 1, 3])
```
"""
function involution(p::SetPartition)
    return set_partition(lower_points(p), upper_points(p))
end

"""
    reflect_vertical(p::SetPartition)

Reflect `p` at the vertical axis.

The vertical reflection of a partition is obtained by reversing the order of 
the upper and lower points. See also Section 4.1.2 in [Gro20](@cite).

# Examples
```jldoctest
julia> reflect_vertical(set_partition([1, 2, 3], [2, 1]))
SetPartition([1, 2, 3], [3, 2])
```
"""
function reflect_vertical(p::SetPartition)
    return set_partition(reverse(upper_points(p)), reverse(lower_points(p)))
end

"""
    rotate_top_left(p::SetPartition)

Perform a top-left rotation of `p`.

A top-left rotation of a partition moves the leftmost point of the upper points 
to the lower points. See also Section 4.1.2 in [Gro20](@cite).

# Examples
```jldoctest
julia> rotate_top_left(set_partition([1, 2, 3], [2, 1]))
SetPartition([1, 2], [3, 1, 3])
```
"""
function rotate_top_left(p::SetPartition)
    
    @req !isempty(upper_points(p)) "partition has no top part"

    result = deepcopy((upper_points(p), lower_points(p)))
    a = result[1][1]
    splice!(result[1], 1)
    pushfirst!(result[2], a)
    
    return set_partition(result[1], result[2])
end

"""
    rotate_bottom_left(p::SetPartition)

Perform a bottom-left rotation of `p`.

A bottom-left rotation of a partition moves the leftmost point of the lower points 
to the upper points. See also Section 4.1.2 in [Gro20](@cite).

# Examples
```jldoctest
julia> rotate_bottom_left(set_partition([1, 2, 3], [2, 1]))
SetPartition([1, 2, 1, 3], [2])
```
"""
function rotate_bottom_left(p::SetPartition)
    
    @req !isempty(lower_points(p)) "partition has no bottom part"

    result = deepcopy((upper_points(p), lower_points(p)))
    a = result[2][1]
    splice!(result[2], 1)
    pushfirst!(result[1], a)
    
    return set_partition(result[1], result[2])
end

"""
    rotate_top_right(p::SetPartition)

Perform a top-right rotation of `p`.

A top-right rotation of a partition moves the rightmost point of the upper points 
to the lower points. See also Section 4.1.2 in [Gro20](@cite).

# Examples
```jldoctest
julia> rotate_top_right(set_partition([1, 2, 3], [2, 1]))
SetPartition([1, 2], [2, 1, 3])
```
"""
function rotate_top_right(p::SetPartition)
    

    @req !isempty(upper_points(p)) "partition has no top part"

    result = deepcopy((upper_points(p), lower_points(p)))
    a = result[1][end]
    pop!(result[1])
    push!(result[2], a)
    
    return set_partition(result[1], result[2])
end

"""
    rotate_bottom_right(p::SetPartition)

Perform a bottom-right rotation of `p`.

A bottom-right rotation of a partition moves the rightmost point of the lower points 
to the upper points. See also Section 4.1.2 in [Gro20](@cite).

# Examples
```jldoctest
julia> rotate_bottom_right(set_partition([1, 2, 3], [2, 1]))
SetPartition([1, 2, 3, 1], [2])
```
"""
function rotate_bottom_right(p::SetPartition)
    

    @req !isempty(lower_points(p)) "partition has no bottom part"

    result = deepcopy((upper_points(p), lower_points(p)))
    a = result[2][end]
    pop!(result[2])
    push!(result[1], a)
    
    return set_partition(result[1], result[2])
end

"""
    is_composable(p::SetPartition, q::SetPartition)

Return whether `p` and `q` are composable, i.e. the number of upper points of 
`p` equals the number of lower points of `q`.

# Examples
```jldoctest
julia> is_composable(set_partition([1, 2], [2, 1]), set_partition([1], [1, 1]))
true

julia> is_composable(set_partition([1], [1, 1]), set_partition([1, 2], [2, 1]))
false
```
"""
function is_composable(p::SetPartition, q::SetPartition)
    return number_of_upper_points(p) == number_of_lower_points(q)
end

"""
    compose_count_loops(p::SetPartition, q::SetPartition)

Return the composition of `p` and `q` as well as the number of removed loops.

The composition of two partitions is obtained by concatenating them vertically
and removing intermediate loops which are no longer connected to the top or bottom.
See also Section 4.1.1 in [Gro20](@cite).

The composition of `p` and `q` is only defined if the number of upper points of 
`p` equals the number of lower points of `q`. See also [`is_composable`](@ref).

# Examples
```jldoctest
julia> compose_count_loops(set_partition([1, 2], [2, 1]), set_partition([1], [1, 1]))
(SetPartition([1], [1, 1]), 0)

julia> compose_count_loops(set_partition([1, 1], [2]), set_partition([1], [2, 2]))
(SetPartition([1], [2]), 1)

julia> compose_count_loops(set_partition([1], [1, 2]), set_partition([1], [2, 2]))
ERROR: ArgumentError: number of points mismatch
[...]
```
"""
function compose_count_loops(p::SetPartition, q::SetPartition)
    
    @req is_composable(p, q) "number of points mismatch" 

    # Work with copies to not change the input partitions
    p_copy = deepcopy(p)

    # new_ids dictionary stores the new value we need to assign to the partition,
    # in order to connect new segments
    vector_q = _new_point_values(upper_points(p_copy), lower_points(p_copy), 
                deepcopy(upper_points(q)), deepcopy(lower_points(q)))
    new_ids = Dict{Int, Int}()
    
    # mapping the second the lower points of the second partition 
    # to the upper points of the first partition and merge if connection
    for (i, n) in enumerate(vector_q[2])
        if !(n in keys(new_ids))
            new_ids[n] = upper_points(p)[i]
        else
            if upper_points(p)[i] in keys(new_ids) && new_ids[n] in keys(new_ids)
                # Do path compression if we have the case that we need to merge two tree's 
                # together and the nodes we operate on are not a root or a leaf
                for ii in [n]
                    path = [ii]
                    z = new_ids[ii]
                    already_in = Set(z)
                    while z in keys(new_ids)
                        push!(path, z)
                        push!(already_in, z)
                        z = new_ids[z]
                        z in already_in && break
                    end
                    push!(path, z)
                    for nn in path[1:end-1]
                        new_ids[nn] = path[end]
                    end
                end
                new_ids[new_ids[n]] = new_ids[upper_points(p_copy)[i]]
            else
                if !(new_ids[n] in keys(new_ids))
                    new_ids[new_ids[n]] = upper_points(p_copy)[i]
                else
                    new_ids[upper_points(p_copy)[i]] = new_ids[n]
                end
            end
        end
    end
    
    # final path compression
    for (ii, z) in new_ids
        path = [ii]
        already_in = Set(z)
        while z in keys(new_ids)
            push!(path, z)
            push!(already_in, z)
            z = new_ids[z]
            z in already_in && break
        end
        push!(path, z)
        for nn in path[1:end-1]
            new_ids[nn] = path[end]
        end
    end
    
    # giving the top part new values
    for (i, n) in enumerate(vector_q[1])
        if n in keys(new_ids)
            vector_q[1][i] = new_ids[n]
        end
    end

    # giving the top part new values
    for (i, n) in enumerate(lower_points(p_copy))
        if n in keys(new_ids)
            lower_points(p_copy)[i] = new_ids[n]
        end
    end

    # removing the middle by just changing the top of our partition 
    # to the adjusted top of the second partition
    ret = set_partition(vector_q[1], lower_points(p_copy))

    # calculating removed related components (loop)
        
    related_comp = Set()
    return_partition_as_set = Set(vcat(vector_q[1], lower_points(p_copy)))

    # calculate new ids for middle nodes, which are under normal circumstances omitted
    for (i, n) in enumerate(vector_q[2])
        if n in keys(new_ids)
            vector_q[2][i] = new_ids[n]
        end
    end

    for (i, n) in enumerate(upper_points(p_copy))
        if n in keys(new_ids)
            upper_points(p_copy)[i] = new_ids[n]
        end
    end

    for co in vcat(vector_q[2], upper_points(p_copy))
        if !(co in return_partition_as_set)
            push!(related_comp, co)
        end
    end
    
    return (ret, length(related_comp))
end

"""
    number_of_blocks(p::SetPartition)

Return the number of blocks of `p`.

# Examples
```jldoctest
julia> number_of_blocks(set_partition([1, 2, 3], [2, 1, 3, 3]))
3
```
"""
function number_of_blocks(p::SetPartition)
    # obtain one vector describing the partition V
    vec = vcat(upper_points(p), lower_points(p))

    # return the maximum number in vec
    maximum(vec; init=0)
end

"""
    is_dominated_by(p::SetPartition, q::SetPartition)

Check if the set partition `p` is dominated by the set partition `q`. This is the case if every block of `p`
is contained in exactly one block of `q`.

# Examples
```jldoctest
julia> is_dominated_by(set_partition([1, 1, 2], [1, 2, 3]), set_partition([1, 1, 2], [1, 2, 1]))
true
```
"""
function is_dominated_by(p::SetPartition, q::SetPartition)
    @req size(p) == size(q) "arguments must have the same size"

    # obtain vectors describing the partitions p and q
    p_vec = vcat(upper_points(p), lower_points(p))
    q_vec = vcat(upper_points(q), lower_points(q))

    # introduce a dictionary to store a mapping from the blocks of p to the blocks of q
    block_map = Dict{Int, Int}()

    for (index, block) in enumerate(p_vec)
        # if the block of the index in p has already been mapped to a block of q,
        # check if the mapping is consistent, otherwise add the mapping
        if haskey(block_map, block) && q_vec[index] != block_map[block]
            return false
        else
            block_map[block] = q_vec[index]
        end
    end
    return true
end

"""
    cycle_partition(p::PermGroupElem)

Return the set partition whose blocks are the cycles of the permutation `p`. This set partition has no lower points.

# Examples
```jldoctest
julia> cycle_partition(perm(symmetric_group(3), [2, 1, 3]))
SetPartition([1, 1, 2], Int64[])
```
"""
function cycle_partition(p::PermGroupElem)
    cycle_list = collect(cycles(p))
    n = degree(parent(p))
    partition_vector = zeros(Int64, n)

    for (index, cycle) in enumerate(cycle_list)
        for element in cycle
            partition_vector[element] = index
        end
    end

    return set_partition(partition_vector, Int64[])
end

"""
    join(p::SetPartition, q::SetPartition)

Return the join of `p` and `q`. This is the unique set partition, where two elements are in the same block of the partition
iff this is the case in `p` or `q`.

# Examples
```jldoctest
julia> join(set_partition([1, 2], [2, 3]), set_partition([1, 2], [1, 3]))
SetPartition([1, 1], [1, 2])
```
"""
function join(p::SetPartition, q::SetPartition)
    @req length(upper_points(p)) == length(upper_points(q)) "p and q must have the same number of upper points"
    @req length(lower_points(p)) == length(lower_points(q)) "p and q must have the same number of lower points"

    _p = set_partition(vcat(upper_points(p), lower_points(p)), vcat(upper_points(p), lower_points(p)))
    _q = set_partition(vcat(upper_points(q), lower_points(q)), vcat(upper_points(q), lower_points(q)))

    join_p_q = upper_points(compose(_p, _q))

    return set_partition(join_p_q[1:length(upper_points(p))], join_p_q[(length(upper_points(p))+1):end])
end
