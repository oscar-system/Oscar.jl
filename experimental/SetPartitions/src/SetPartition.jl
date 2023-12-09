"""
    SetPartition

Initialize SetPartition object, while transforming object to "normal form". 

# Examples
```jldoctest
julia> SetPartition([2, 4], [4, 99])
SetPartition([1, 2], [2, 3])
```
"""
struct SetPartition <: AbstractPartition
    upper_points::Vector{Int}
    lower_points::Vector{Int}

    function SetPartition(_upper_points, _lower_points)
        (__upper_points, __lower_points) = 
            normal_form_vector(copy([Int.(copy(_upper_points)), 
                Int.(copy(_lower_points))]))
        new(__upper_points, __lower_points)
    end
end


function hash(p::SetPartition, h::UInt)
    hash(p.upper_points, hash(p.lower_points, h))
end

function ==(p::SetPartition, q::SetPartition)

    p.lower_points == q.lower_points && p.upper_points == q.upper_points

end

function copy(p::SetPartition)
    SetPartition(copy(p.upper_points), copy(p.lower_points))
end

"""
    tensor_product(p::SetPartition, q::SetPartition)

Return tensor product of `p` and `q`.

# Examples
```jldoctest
julia> tensor_product(SetPartition([1, 2], [2, 1]), SetPartition([1, 1], [1]))
SetPartition([1, 2, 3, 3], [2, 1, 3])
```
"""
function tensor_product(p::SetPartition, q::SetPartition)
    
    q_new = helper_new_point_values_vector([p.upper_points, p.lower_points], 
        [q.upper_points, q.lower_points])
    SetPartition(vcat(p.upper_points, q_new[1]), vcat(p.lower_points, q_new[2]))
end

"""
    involution(p::SetPartition)

Return involution of `p`.

# Examples
```jldoctest
julia> involution(SetPartition([1, 2, 3], [2, 1]))
SetPartition([1, 2], [2, 1, 3])
```
"""
function involution(p::SetPartition)
    
    SetPartition(p.lower_points, p.upper_points)

end

"""
    vertical_reflection(p::SetPartition)

Return vertical reflection of `p`.

# Examples
```jldoctest
julia> vertical_reflection(SetPartition([1, 2, 3], [2, 1]))
SetPartition([1, 2, 3], [3, 2])
```
"""
function vertical_reflection(p::SetPartition)
    
    SetPartition(reverse(p.upper_points), reverse(p.lower_points))

end

"""
    rotation(p::SetPartition)

Return the rotation of `p` in the direction given by `lr` and `tb`. 

# Arguments
- `p`: input partition
- `lr`: whether left (true) or right (false)
- `tb`: whether top (true) or bottom (false)

# Examples
```jldoctest
julia> rotation(SetPartition([1, 2, 3], [2, 1]), true, true)
SetPartition([1, 2], [3, 1, 3])
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
    
    SetPartition(ret[1], ret[2])
end

"""
    composition_loops(p::SetPartition, q::SetPartition)

Return the composition of `p` and `q` as well as the number of removed loops.

# Examples
```jldoctest
julia> composition_loops(SetPartition([1, 2], [2, 1]), SetPartition([1], [1, 1]))
(SetPartition([1], [1, 1]), 0)
```
"""
function composition_loops(p::SetPartition, q::SetPartition)
    
    !is_composable(p, q) ? error("format not fitting") : 

    # Work with copies to not change the input partitions
    p_copy = copy(p)

    # new_ids dictionary stores the new value we need to assign to the partition,
    # in order to connect new segments
    vector_q = helper_new_point_values_vector([p_copy.upper_points, p_copy.lower_points], 
        [copy(q.upper_points), copy(q.lower_points)])
    new_ids = Dict{Int, Int}()
    
    # mapping the second the lower points of the second partition 
    # to the upper points of the first partition and merge if connection
    for (i, n) in enumerate(vector_q[2])
        if !(n in keys(new_ids))
            new_ids[n] = p.upper_points[i]
        else
            if p.upper_points[i] in keys(new_ids) && get(new_ids, n, -1) in keys(new_ids)
                # Do path compression if we have the case that we need to merge two tree's 
                # together and the nodes we operate on are not a root or a leaf
                for ii in [n]
                    path = [ii]
                    z = get(new_ids, ii, -1)
                    already_in = Set(z)
                    while z in keys(new_ids)
                        push!(path, z)
                        push!(already_in, z)
                        z = get(new_ids, z, -1)
                        z in already_in && break
                    end
                    push!(path, z)
                    for nn in path[1:end-1]
                        new_ids[nn] = path[end]
                    end
                end
                new_ids[get(new_ids, n, -1)] = get(new_ids, p_copy.upper_points[i], -1)
            else
                if !(get(new_ids, n, -1) in keys(new_ids))
                    new_ids[get(new_ids, n, -1)] = p_copy.upper_points[i]
                else
                    new_ids[p_copy.upper_points[i]] = get(new_ids, n, 1)
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
            z = get(new_ids, z, -1)
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
            vector_q[1][i] = get(new_ids, n, -1)
        end
    end

    # giving the top part new values
    for (i, n) in enumerate(p_copy.lower_points)
        if n in keys(new_ids)
            p_copy.lower_points[i] = get(new_ids, n, -1)
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
            vector_q[2][i] = get(new_ids, n, -1)
        end
    end

    for (i, n) in enumerate(p_copy.upper_points)
        if n in keys(new_ids)
            p_copy.upper_points[i] = get(new_ids, n, -1)
        end
    end

    for co in vcat(vector_q[2], p_copy.upper_points)
        if !(co in return_partition_as_set)
            push!(related_comp, co)
        end
    end
    
    (ret, length(related_comp))
end
