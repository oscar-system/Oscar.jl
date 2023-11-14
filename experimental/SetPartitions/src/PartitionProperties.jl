function size(p::SetPartition)
    length(p.lower_points) + length(p.upper_points)
end

function is_composable(p::SetPartition, q::SetPartition)
    length(p.upper_points) == length(q.lower_points)
end

function upper_points(p::SetPartition)
    p.upper_points
end

function lower_points(p::SetPartition)
    p.lower_points
end

function size(p::ColoredPartition)
    size(p.partition)
end

function is_composable(p::ColoredPartition, q::ColoredPartition)
    p.color_upper_points == q.color_lower_points && is_composable(p.partition, q.partition)
end

function upper_points(p::ColoredPartition)
    p.partition.upper_points
end

function lower_points(p::ColoredPartition)
    p.partition.lower_points
end

function size(p::SpatialPartition)
    size(p.partition)
end

function is_composable(p::SpatialPartition, q::SpatialPartition)
    p.dimension == q.dimension && is_composable(p.partition, q.partition)
end

function is_worth_composition(p::AbstractPartition, q::AbstractPartition, max_length::Int)
    length(upper_points(p)) != 0 && length(upper_points(p)) != max_length && 
        length(lower_points(p)) + length(upper_points(q)) <= max_length
end

function upper_points(p::SpatialPartition)
    p.partition.upper_points
end

function lower_points(p::SpatialPartition)
    p.partition.lower_points
end

"""
    is_pair(p::AbstractPartition)

Check and output whether `p` is a partition only consisting of blocks of size two.

# Examples
```jldoctest
julia> is_pair(Partition([1, 2, 2, 1, 3], [3]))
true
```
"""
function is_pair(p::AbstractPartition)

    # Dictionary from block to size of block
    block_to_size = Dict{Int, Int}()

    # Initialize dictionary
    for i in vcat(upper_points(p), lower_points(p))
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
    is_balanced(p::T) where {T<:Union{SetPartition, ColoredPartition}}

Check and output whether `p` is a balanced partition.

# Examples
```jldoctest
julia> is_balanced(SetPartition([1, 2, 3], [3, 2, 1]))
true
```
"""
function is_balanced(p::T) where {T<:Union{SetPartition, ColoredPartition}}

    p_vector = vcat(upper_points(p), lower_points(p))

    # remember length lower points
    upper = length(upper_points(p))
    upper_uneven = upper % 2 == 1 ? true : false

    # Dictionary from block to sum of -1 (repr odd indices) and 1 (repr even indices)
    block_to_size = Dict{Int, Int}()

    # prefill dict with zeros
    for i in p_vector
        block_to_size[i] = 0
    end

    # Initialize dictionary
    for (i, n) in enumerate(p_vector)
        condition = upper_uneven ? i % 2 == 1 : (i <= upper ? i % 2 == 1 : i % 2 == 0)
        if condition
            block_to_size[n] = get(block_to_size, n, -1) - 1
        else
            block_to_size[n] = get(block_to_size, n, -1) + 1
        end
    end

    return all(i -> i == 0, values(block_to_size))
end

"""
    is_noncrossing(p::T) where {T<:Union{SetPartition, ColoredPartition}}

Check and output whether `p` is a non-crossing partition.

# Examples
```jldoctest
julia> is_noncrossing(SetPartition([1, 2, 2, 3, 1, 4], [4, 3]))
false
```
"""
function is_noncrossing(p::T) where {T<:Union{SetPartition, ColoredPartition}}

    # transform partition to only upper points
    p_vector = vcat(upper_points(p), reverse(lower_points(p)))

    # Dictionary from block to size of block
    block_to_size = Dict{Int, Int}()

    # Initialize dictionary
    for i in p_vector
        if !(i in keys(block_to_size))
            block_to_size[i] = 1
        else
            block_to_size[i] = get(block_to_size, i, -1) + 1
        end
    end

    # blocks we have already seen in the iteration process
    already_seen = Set{Int}()
    last_incompleted = []
    incompleted = Set{Int}()

    for (i, n) in enumerate(p_vector)
        if (n in incompleted && 
            (isempty(last_incompleted) ? -1 : last_incompleted[end]) != n)
            return false
        end
        if !(n in already_seen)
            push!(already_seen, n)
            block_to_size[n] = get(block_to_size, n, -1) - 1
            if (get(block_to_size, n, -1) > 0 && 
                (isempty(last_incompleted) ? -1 : last_incompleted[end]) != n)
                push!(last_incompleted, n)
                push!(incompleted, n)
            end
        else
            if p_vector[i-1] != n && get(block_to_size, p_vector[i-1], -1) != 0
                return false
            else
                block_to_size[n] = get(block_to_size, n, -1) - 1
                if get(block_to_size, n, -1) == 0 && !isempty(last_incompleted)
                    if last_incompleted[end] == n
                        pop!(last_incompleted)
                    end
                    delete!(incompleted, n)
                end
                if get(block_to_size, n, -1) > 0 && 
                    (isempty(last_incompleted) ? -1 : last_incompleted[end]) != n
                    push!(last_incompleted, n)
                    push!(incompleted, n)
                end
            end
        end
    end
    return true
end
