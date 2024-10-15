"""
    is_pair(p::AbstractPartition)

Return whether `p` is a partition consisting only of subsets of size two.

# Examples
```jldoctest
julia> is_pair(set_partition([1, 2, 2, 1, 3], [3]))
true

julia> is_pair(set_partition([1, 2, 3, 1, 3], [3]))
false
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
            block_to_size[i] = block_to_size[i] + 1
            if block_to_size[i] > 2
                return false
            end
        end
    end
    return all(==(2), values(block_to_size))
end

"""
    is_balanced(p::T) where {T<:Union{SetPartition, ColoredPartition}}

Return whether `p` is a balanced partition, i.e. every subset has the 
same number of points with even and odd indices.

Note that the parity of the lower indices is inverted.

# Examples
```jldoctest
julia> is_balanced(set_partition([1, 2, 3], [3, 2, 1]))
true

julia> is_balanced(set_partition([1, 2, 3], [3, 1, 2]))
false
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
            block_to_size[n] = block_to_size[n] - 1
        else
            block_to_size[n] = block_to_size[n] + 1
        end
    end

    return all(iszero, values(block_to_size))
end

"""
    is_non_crossing(p::T) where {T<:Union{SetPartition, ColoredPartition}}

Return whether `p` is a non-crossing partition, i.e. there are no crossing
lines when drawing `p`. See also Section 4.1.1 in [Gro20](@cite).

# Examples
```jldoctest
julia> is_non_crossing(set_partition([1, 2, 2, 1, 3], [3, 3]))
true

julia> is_non_crossing(set_partition([1, 2, 2, 3, 1, 4], [4, 3]))
false
```
"""
function is_non_crossing(p::T) where {T<:Union{SetPartition, ColoredPartition}}

    # transform partition to only upper points
    p_vector = vcat(upper_points(p), reverse(lower_points(p)))

    # Dictionary from block to size of block
    block_to_size = Dict{Int, Int}()

    # Initialize dictionary
    for i in p_vector
        block_to_size[i] = get(block_to_size, i, 0) + 1
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
            block_size = block_to_size[n] = block_to_size[n] - 1
            if (block_size > 0 && 
                (isempty(last_incompleted) ? -1 : last_incompleted[end]) != n)
                push!(last_incompleted, n)
                push!(incompleted, n)
            end
        else
            if p_vector[i-1] != n && get(block_to_size, p_vector[i-1], -1) != 0
                return false
            end
            block_size = block_to_size[n] = block_to_size[n] - 1
            if block_size == 0 && !isempty(last_incompleted)
                if last_incompleted[end] == n
                    pop!(last_incompleted)
                end
                delete!(incompleted, n)
            end
            if block_size > 0 && 
                (isempty(last_incompleted) ? -1 : last_incompleted[end]) != n
                push!(last_incompleted, n)
                push!(incompleted, n)
            end
        end
    end
    return true
end
