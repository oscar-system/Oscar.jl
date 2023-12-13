"""
    SetPartition

A set-partition is a partition of a set of upper and lower points into disjoint subsets. 
Such partitions are often depicted as string diagrams where points in the same subset 
are connected by lines. See also Section 4.1.1 in [Gro20](@cite).

`SetPartition` stores two vectors `upper_points` and `lower_points`
where same numbers correspond to points in the same subset. Both vectors are 
automatically converted into normal form as returned by [`normal_form`](@ref).

# Examples
```jldoctest
julia> SetPartition([2, 4], [4, 99])
SetPartition([1, 2], [2, 3])
```
"""
struct SetPartition <: AbstractPartition
    upper_points::Vector{Int}
    lower_points::Vector{Int}

    function SetPartition(upper_points, lower_points)
        (new_upper, new_lower) = normal_form(convert(Vector{Int}, upper_points), 
                                             convert(Vector{Int}, lower_points))
        return new(new_upper, new_lower)
    end
end


function hash(p::SetPartition, h::UInt)
    return hash(p.upper_points, hash(p.lower_points, h))
end

function ==(p::SetPartition, q::SetPartition)
    return p.lower_points == q.lower_points && p.upper_points == q.upper_points
end

function copy(p::SetPartition)
    return SetPartition(copy(p.upper_points), copy(p.lower_points))
end

"""
    tensor_product(p::SetPartition, q::SetPartition)

Return tensor product of `p` and `q`.

The tensor product of two partitions is given by their vertical concatenation.
See also Section 4.1.1 in [Gro20](@cite).

# Examples
```jldoctest
julia> tensor_product(SetPartition([1, 2], [2, 1]), SetPartition([1, 1], [1]))
SetPartition([1, 2, 3, 3], [2, 1, 3])
```
"""
function tensor_product(p::SetPartition, q::SetPartition)
    
    q_new = new_point_values(p.upper_points, p.lower_points, 
                             q.upper_points, q.lower_points)
    return SetPartition(vcat(p.upper_points, q_new[1]), vcat(p.lower_points, q_new[2]))
end

"""
    involution(p::SetPartition)

Return involution of `p`.

The involution of a partition is obtained by swapping the upper and lower points.
See also Section 4.1.1 in [Gro20](@cite).

# Examples
```jldoctest
julia> involution(SetPartition([1, 2, 3], [2, 1]))
SetPartition([1, 2], [2, 1, 3])
```
"""
function involution(p::SetPartition)
    return SetPartition(p.lower_points, p.upper_points)
end

"""
    vertical_reflection(p::SetPartition)

Return vertical reflection of `p`.

The vertical reflection of a partition is obtained by reversing the order of 
the upper and lower points. See also Section 4.1.2 in [Gro20](@cite).

# Examples
```jldoctest
julia> vertical_reflection(SetPartition([1, 2, 3], [2, 1]))
SetPartition([1, 2, 3], [3, 2])
```
"""
function vertical_reflection(p::SetPartition)
    return SetPartition(reverse(p.upper_points), reverse(p.lower_points))
end

"""
    rotation(p::SetPartition)

Return the rotation of `p` in the direction given by `lr` and `tb`. 

Rotating a partition moves the left- or right-most point of the upper points 
to the lower points or vice verca. See also Section 4.1.2 in [Gro20](@cite).

# Arguments
- `p`: input partition
- `lr`: rotating at the left (true) or at the right (false)
- `tb`: rotating from top to bottom (true) or from bottom to top (false)

# Examples
```jldoctest
julia> rotation(SetPartition([1, 2, 3], [2, 1]), true, true)
SetPartition([1, 2], [3, 1, 3])

julia> rotation(SetPartition([1, 2, 3], [2, 1]), true, false)
SetPartition([1, 2, 1, 3], [2])

julia> rotation(SetPartition([1, 2, 3], [2, 1]), false, true)
SetPartition([1, 2], [2, 1, 3])

julia> rotation(SetPartition([1, 2, 3], [2, 1]), false, false)
SetPartition([1, 2, 3, 1], [2])
```
"""
function rotation(p::SetPartition, lr::Bool, tb::Bool)
    
    if tb && isempty(p.upper_points)
        error("SetPartition has no top part")
    elseif !tb && isempty(p.lower_points)
        error("SetPartition has no bottom part")
    end

    ret = (copy(p.upper_points), copy(p.lower_points))

    if lr
        if tb
            a = ret[1][1]
            splice!(ret[1], 1)
            pushfirst!(ret[2], a)
        else
            a = ret[2][1]
            splice!(ret[2], 1)
            pushfirst!(ret[1], a)
        end
    else
        if tb
            a = ret[1][end]
            pop!(ret[1])
            push!(ret[2], a)
        else
            a = ret[2][end]
            pop!(ret[2])
            push!(ret[1], a)
        end
    end
    
    return SetPartition(ret[1], ret[2])
end

"""
    composition_loops(p::SetPartition, q::SetPartition)

Return the composition of `p` and `q` as well as the number of removed loops.

The composition of two partitions is obtained concatenating them horizontally
and removing intermediate loops which are no longer connected to the top or bottom.
See also Section 4.1.1 in [Gro20](@cite).

The composition of `p` and `q` is only defined if the number of upper points of 
`p` equals the number of lower points of `q`. See also [`is_composable`](@ref).

# Examples
```jldoctest
julia> composition_loops(SetPartition([1, 2], [2, 1]), SetPartition([1], [1, 1]))
(SetPartition([1], [1, 1]), 0)

julia> composition_loops(SetPartition([1, 1], [2]), SetPartition([1], [2, 2]))
(SetPartition([1], [2]), 1)

julia> composition_loops(SetPartition([1], [1, 2]), SetPartition([1], [2, 2]))
ERROR: format not fitting
```
"""
function composition_loops(p::SetPartition, q::SetPartition)
    
    !is_composable(p, q) ? error("format not fitting") : 

    # Work with copies to not change the input partitions
    p_copy = copy(p)

    # new_ids dictionary stores the new value we need to assign to the partition,
    # in order to connect new segments
    vector_q = new_point_values(p_copy.upper_points, p_copy.lower_points, 
                                copy(q.upper_points), copy(q.lower_points))
    new_ids = Dict{Int, Int}()
    
    # mapping the second the lower points of the second partition 
    # to the upper points of the first partition and merge if connection
    for (i, n) in enumerate(vector_q[2])
        if !(n in keys(new_ids))
            new_ids[n] = p.upper_points[i]
        else
            if p.upper_points[i] in keys(new_ids) && new_ids[n] in keys(new_ids)
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
                new_ids[new_ids[n]] = new_ids[p_copy.upper_points[i]]
            else
                if !(new_ids[n] in keys(new_ids))
                    new_ids[new_ids[n]] = p_copy.upper_points[i]
                else
                    new_ids[p_copy.upper_points[i]] = new_ids[n]
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
    for (i, n) in enumerate(p_copy.lower_points)
        if n in keys(new_ids)
            p_copy.lower_points[i] = new_ids[n]
        end
    end

    # removing the middle by just changing the top of our partition 
    # to the adjusted top of the second partition
    ret = SetPartition(vector_q[1], p_copy.lower_points)

    # calculating removed related components (loop)
        
    related_comp = Set()
    return_partition_as_set = Set(vcat(vector_q[1], p_copy.lower_points))

    # calculate new ids for middle nodes, which are under normal circumstances omitted
    for (i, n) in enumerate(vector_q[2])
        if n in keys(new_ids)
            vector_q[2][i] = new_ids[n]
        end
    end

    for (i, n) in enumerate(p_copy.upper_points)
        if n in keys(new_ids)
            p_copy.upper_points[i] = new_ids[n]
        end
    end

    for co in vcat(vector_q[2], p_copy.upper_points)
        if !(co in return_partition_as_set)
            push!(related_comp, co)
        end
    end
    
    return (ret, length(related_comp))
end
