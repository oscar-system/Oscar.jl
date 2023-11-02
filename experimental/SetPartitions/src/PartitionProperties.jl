"""
size(p::Partition)

This function outputs the size of the input partition (i.e. the number of points in `p`).

# Arguments
- `p`: Input partition

# Returns
- size of `p`
"""
function size(p::Partition)
    length(p.lower_points) + length(p.upper_points)
end

function is_composable(p::Partition, q::Partition)
    length(p.upper_points) == length(q.lower_points)
end

function get_upper_points(p::Partition)
    p.upper_points
end

function get_lower_points(p::Partition)
    p.lower_points
end

"""
size(p::ColoredPartition)

This function outputs the size of the input colroed partition (i.e. the number of points in `p`).

# Arguments
- `p`: Input colored partition

# Returns
- size of `p`
"""
function size(p::ColoredPartition)
    size(p.partition)
end

function is_composable(p::ColoredPartition, q::ColoredPartition)
    p.color_upper_points == q.color_lower_points && is_composable(p.partition, q.partition)
end

function get_upper_points(p::ColoredPartition)
    p.partition.upper_points
end

function get_lower_points(p::ColoredPartition)
    p.partition.lower_points
end

"""
size(p::SpatialPartition)

This function outputs the size of the input spatial partition (i.e. the number of points in `p`).

# Arguments
- `p`: Input spatial partition

# Returns
- size of `p`
"""
function size(p::SpatialPartition)
    size(p.partition)
end

function is_composable(p::SpatialPartition, q::SpatialPartition)
    p.dimension == q.dimension && is_composable(p.partition, q.partition)
end

function is_worth_composition(p::AbstractPartition, q::AbstractPartition, max_length::Int)
    length(get_upper_points(p)) != 0 && length(get_upper_points(p)) != max_length && length(get_lower_points(p)) + length(get_upper_points(q)) <= max_length
end

function get_upper_points(p::SpatialPartition)
    p.partition.upper_points
end

function get_lower_points(p::SpatialPartition)
    p.partition.lower_points
end

"""
is_pair(p::Partition)

This function checks whether `p` is a partition only including blocks of size two (in O(n) average).

# Arguments
- `p`: Input partition

# Returns
- true if `p` is pair partition else false

# Examples
```julia-repl
julia> is_pair(Partition([1, 2, 2, 1, 3], [3]))
true
```
"""
function is_pair(p::AbstractPartition)

    # Dictionary from block to size of block
    block_to_size = Dict()

    # Initialize dictionary
    for i in vcat(get_upper_points(p), get_lower_points(p))
        if !(i in keys(block_to_size))
            block_to_size[i] = 1
        else
            block_to_size[i] = get(block_to_size, i, -1) + 1
            if get(block_to_size, i, -1) > 2
                return false
            end
        end
    end
    return all(i -> i == 2, values(block_to_size))
end

"""
is_balanced(p::Partition)

This function checks whether `p` is a balanced partition (in O(n) average).

# Arguments
- `p`: Input partition

# Returns
- true if `p` is balanced partition else false

# Examples
```julia-repl
julia> is_balanced(Partition([1, 2, 3], [3, 2, 1]))
true
```
"""
function is_balanced(p::T) where {T<:Union{Partition, ColoredPartition}}

    p_array = vcat(get_upper_points(p), get_lower_points(p))

    # remember length lower points
    upper = length(get_upper_points(p))
    upper_uneven = upper % 2 == 1 ? true : false

    # Dictionary from block to sum of -1 (repr odd indices) and 1 (repr even indices)
    block_to_size = Dict()
    block_to_blocksize = Dict()

    # prefill dict with zeros
    for i in p_array
        block_to_size[i] = 0
    end
    # to identify singeltons
    # for i in p_array
    #     if i in keys(block_to_blocksize)
    #        block_to_blocksize[i] = get(block_to_blocksize, i, -1) + 1
    #    else
    #        block_to_blocksize[i] = 1
    #    end
    #end

    # Initialize dictionary
    for (i, n) in enumerate(p_array)
        condition = upper_uneven ? i % 2 == 1 : (i <= upper ? i % 2 == 1 : i % 2 == 0)
        if condition # && get(block_to_blocksize, n, -1) != 1
            block_to_size[n] = get(block_to_size, n, -1) - 1
        else # if get(block_to_blocksize, n, -1) != 1
            block_to_size[n] = get(block_to_size, n, -1) + 1
        end
    end
    # check whether singeltons are balanced
    #singeltons = 0
    #for i in keys(block_to_blocksize)
    #    if get(block_to_blocksize, i, -1) != 1
    #        continue
    #    end
    #    i % 2 == 1 ? singeltons -= 1 : singeltons += 1
    #end

    return all(i -> i == 0, values(block_to_size)) #&& singeltons == 0
end

"""
is_noncrossing(p::Partition)

This function checks whether `p` is a non-crossing partition (in O(n) average).

# Arguments
- `p`: Input partition

# Returns
- true if `p` is non-crossing partition else false

# Examples
```julia-repl
julia> is_noncrossing(Partition([1, 2, 2, 3, 1, 4], [4, 3]))
false
```
"""
function is_noncrossing(p::T) where {T<:Union{Partition, ColoredPartition}}

    # transform partition to only upper points
    p_array = vcat(get_upper_points(p), reverse(get_lower_points(p)))

    # Dictionary from block to size of block
    block_to_size = Dict()

    # Initialize dictionary
    for i in p_array
        if !(i in keys(block_to_size))
            block_to_size[i] = 1
        else
            block_to_size[i] = get(block_to_size, i, -1) + 1
        end
    end

    # blocks we have already seen in the iteration process
    already_seen = Set()
    last_incompleted = []
    incompleted = Set()

    for (i, n) in enumerate(p_array)
        if n in incompleted && (isempty(last_incompleted) ? -1 : last_incompleted[end]) != n
            return false
        end
        if !(n in already_seen)
            push!(already_seen, n)
            block_to_size[n] = get(block_to_size, n, -1) - 1
            if get(block_to_size, n, -1) > 0 && (isempty(last_incompleted) ? -1 : last_incompleted[end]) != n
                push!(last_incompleted, n)
                push!(incompleted, n)
            end
        else
            if p_array[i-1] != n && get(block_to_size, p_array[i-1], -1) != 0
                return false
            else
                block_to_size[n] = get(block_to_size, n, -1) - 1
                if get(block_to_size, n, -1) == 0 && !isempty(last_incompleted)
                    if last_incompleted[end] == n
                        pop!(last_incompleted)
                    end
                    delete!(incompleted, n)
                end
                if get(block_to_size, n, -1) > 0 && (isempty(last_incompleted) ? -1 : last_incompleted[end]) != n
                    push!(last_incompleted, n)
                    push!(incompleted, n)
                end
            end
        end
    end
    return true
end