import Base.hash
import Base.==
import Base.copy

"""
Partition

Initialize Partition object, while transforming object to "normal form". 

# Arguments
- `upper_points`: upper points of Partition as array
- `lower_points`: lower points of Partition as array

# Examples
```julia-repl
julia> Partition([2, 4], [4, 99])
Partition([1, 2], [2, 3])
```
"""
struct Partition <: AbstractPartition
    upper_points::Array{Int64, 1}
    lower_points::Array{Int64, 1}

    function Partition(_upper_points, _lower_points)
        (__upper_points, __lower_points) = normal_form_array(copy([Int64.(copy(_upper_points)), Int64.(copy(_lower_points))]))
        a = new(__upper_points, __lower_points)
        return a
    end
end


function hash(p::Partition, h::UInt)
    hash(p.upper_points, hash(p.lower_points, h))
end

function ==(p::Partition, q::Partition)

    p.lower_points == q.lower_points && p.upper_points == q.upper_points

end

function copy(p::Partition)
    return Partition(copy(p.upper_points), copy(p.lower_points))
end

"""
tensor_product(p::Partition, q::Partition)

This function applies on p tensor product with q (in O(n)).

# Arguments
- `p`: Input partition
- `q`: Second input partition

# Returns
- `p` tensor product `q`

# Examples
```julia-repl
julia> tensor_product(Partition([1, 2], [2, 1]), Partition([1, 1], [1]))
Partition([1, 2, 3, 3], [2, 1, 3])
```
"""
function tensor_product(p::Partition, q::Partition)
    
    q_new = helper_new_point_values_array([p.upper_points, p.lower_points], [q.upper_points, q.lower_points])
    Partition(vcat(p.upper_points, q_new[1]), vcat(p.lower_points, q_new[2]))
end

"""
involution(p::Partition)

This function applies an involution on `p` (in O(n) because normal_form, else O(1)).

# Arguments
- `p`: Input partition

# Returns
- involution of `p`

# Examples
```julia-repl
julia> involution(Partition([1, 2, 3], [2, 1]))
Partition([1, 2], [2, 1, 3])
```
"""
function involution(p::Partition)
    
    Partition(copy(p.lower_points), copy(p.upper_points))

end

"""
vertical_reflection(p::Partition)

This function applies an vertical reflection on `p` (in O(n) because normal_form, else O(1)).

# Arguments
- `p`: Input partition

# Returns
- vertical reflection of `p`

# Examples
```julia-repl
julia> vertical_reflection(Partition([1, 2, 3], [2, 1]))
Partition([1, 2, 3], [3, 2])
```
"""
function vertical_reflection(p::Partition)
    
    Partition(reverse(p.upper_points), reverse(p.lower_points))

end

"""
rotation(p::Partition)

This function applies a rotation on `p` (in O(n) because normal_form, else O(1)). 
Throws error if rotation not possible.

# Arguments
- `p`: Input partition
- `lr`: lr whether left (true) or right (false)
- `tb`: tb whether top (true) or bottom (false) rotation

# Returns
- rotation of `p`

# Examples
```julia-repl
julia> rotation(Partition([1, 2, 3], [2, 1]), true, true)
Partition([1, 2], [3, 1, 3])
```
"""
function rotation(p::Partition, lr::Bool, tb::Bool)
    
    if tb && isempty(p.upper_points)
        error("Partition has no top part.")
    elseif !tb && isempty(p.lower_points)
        error("Partition has no bottom part.")
    end

    ret::Array = [copy(p.upper_points), copy(p.lower_points)]

        if lr
            if tb
                a::Int = ret[1][1]
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
    
    Partition(ret[1], ret[2])
end

"""
composition_loops(p::Partition, q::Partition)

This function applies composition between p and q (in O(nlogn)).

# Arguments
- `p`: Input partition
- `q`: Second input partition

# Returns
- [`p` composition `q`, number of loops]

# Examples
```julia-repl
julia> composition_loops(Partition([1, 2], [2, 1]), Partition([1], [1, 1]))
Partition([1], [1, 1]), 0
```
"""
function composition_loops(p::Partition, q::Partition)
    
    !is_composable(p, q) ? error("format not fitting") : 

    # Work with copies to not change the input partitions
    p_copy::Partition = copy(p)

    # new_ids dicts store the new Value we need to assign to the partition in order to connect new segments
    array_q = helper_new_point_values_array([p_copy.upper_points, p_copy.lower_points], [copy(q.upper_points), copy(q.lower_points)])
    new_ids::Dict{Int, Int} = Dict()
    
    # fitting the second partition-values to the first and changing if connection
    for (i, n) in enumerate(array_q[2])
        if !(n in keys(new_ids))
            new_ids[n] = p.upper_points[i]
        else
            if p.upper_points[i] in keys(new_ids) && get(new_ids, n, -1) in keys(new_ids)
                # Do path compression if we have the case that we need to merge two tree's together and
                # the nodes we operate on are not a root or a leaf
                for ii::Int in [n]
                    path::Array = [ii]
                    already_in = Set()
                    push!(already_in, get(new_ids, ii, -1))
                    z::Int = get(new_ids, ii, -1)
                    while z::Int in keys(new_ids)
                        push!(path, z)
                        push!(already_in, z)
                        z = get(new_ids, z, -1)
                        if z in already_in
                            break
                        end
                    end
                    push!(path, z)
                    for nn::Int in path[1:end-1]
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
    for ii in keys(new_ids)
        path = [ii]
        already_in = Set()
        push!(already_in, get(new_ids, ii, -1))
        z = get(new_ids, ii, -1)
        while z in keys(new_ids)
            push!(path, z)
            push!(already_in, z)
            z = get(new_ids, z, -1)
            if z in already_in
                break
            end
        end
        push!(path, z)
        for nn in path[1:end-1]
            new_ids[nn] = path[end]
        end
    end
    
    # giving the top part new values
    for (i, n) in enumerate(array_q[1])
        if n in keys(new_ids)
            array_q[1][i] = get(new_ids, n, -1)
        end
    end

    # giving the top part new values
    for (i, n) in enumerate(p_copy.lower_points)
        if n in keys(new_ids)
            p_copy.lower_points[i] = get(new_ids, n, -1)
        end
    end

    # removing the middle by just changing the top of our partition to the adjusted top of the second partition
    ret = Partition(array_q[1], p_copy.lower_points)

    # calculating removed related components (loop)
        
    related_comp = Set()
    return_partition_as_set = Set(vcat(array_q[1], p_copy.lower_points))

    # calculate new ids for middle nodes, which are under normal circumstances omitted
    for (i, n) in enumerate(array_q[2])
        if n in keys(new_ids)
            array_q[2][i] = get(new_ids, n, -1)
        end
    end

    for (i, n) in enumerate(p_copy.upper_points)
        if n in keys(new_ids)
            p_copy.upper_points[i] = get(new_ids, n, -1)
        end
    end

    # if there is a ID in the middle part which is not in result partition set we know, that this is a loop
    for co in vcat(array_q[2], p_copy.upper_points)
        if !(co in return_partition_as_set)
            push!(related_comp, co)
        end
    end
    
    return [ret, length(related_comp)]
end