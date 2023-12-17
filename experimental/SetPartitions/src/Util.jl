"""
    _new_point_values(p_upper::Vector{Int}, p_lower::Vector{Int}, 
                      q_upper::Vector{Int}, q_lower::Vector{Int})

Return two vectors which are semantically identical to `q_upper` and `q_lower`
as set-partition and have no common values with `p_upper` and `p_lower`.
"""
function _new_point_values(p_upper::Vector{Int}, p_lower::Vector{Int}, 
                          q_upper::Vector{Int}, q_lower::Vector{Int})

    p_points = vcat(p_upper, p_lower)

    if !(isempty(p_points))
        new_id, index = findmax(p_points)
        new_id += 1
    else
        new_id = 1 
    end
    
    q_points = vcat(q_upper, q_lower)
    new_ids = Dict{Int, Int}()

    for n in q_points
        new_ids[n] = new_id
        new_id += 1
    end
    
    upper = [new_ids[n] for n in q_upper]
    lower = [new_ids[n] for n in q_lower]
    
    return (upper, lower)
end

"""
    _normal_form(p_upper::Vector{Int}, p_lower::Vector{Int})

Return two vectors which are semantically identical to `p_upper` and `p_lower`
as set-partition and are numbered from 1 to the number of subsets.
"""
function _normal_form(p_upper::Vector{Int}, p_lower::Vector{Int})

    new_id = 1
    new_ids = Dict{Int, Int}()
    p_return = _new_point_values([length(p_upper) + length(p_lower)], 
        Vector{Int}(), p_upper, p_lower)

    for (i, n) in enumerate(p_return[1])
        if !(n in keys(new_ids))
            new_ids[n] = new_id
            p_return[1][i] = new_id
            new_id += 1
        else
            p_return[1][i] = new_ids[n]
        end
    end

    for (i, n) in enumerate(p_return[2])
        if !(n in keys(new_ids))
            new_ids[n] = new_id
            p_return[2][i] = new_id
            new_id += 1
        else
            p_return[2][i] = new_ids[n]
        end
    end

    return p_return
end

"""
    _req_colored_partition(partition::AbstractPartition, 
                            color_upper_points::Vector{Int}, 
                            color_lower_points::Vector{Int})

Return true if the format of the attributes for initializing a
colored partition is valid, false otherwise.

Note that this is a helper function for the constructor of colored partitions.
"""
function _req_colored_partition(partition::AbstractPartition, 
                                color_upper_points::Vector{Int}, 
                                color_lower_points::Vector{Int})
    return all(x -> x in (0, 1), Set(color_upper_points)) &&
            all(x -> x in (0, 1), Set(color_lower_points)) &&
            length(upper_points(partition)) == length(color_upper_points) && 
            length(lower_points(partition)) == length(color_lower_points)
end

"""
    _req_spatial_partition(partition::AbstractPartition, 
                            levels::Int)

Return true if the format of the attributes for initializing a
spatial partition is valid, false otherwise.

Note that this is a helper function for the constructor of spatial partitions.
"""
function _req_spatial_partition(partition::AbstractPartition, 
                                    levels::Int)
    return length(upper_points(partition)) % levels == 0 &&
            length(lower_points(partition)) % levels == 0 &&
            levels > 0
end

"""
    _is_worth_composition(p::T, q::T, max_length::Int) where {T <: AbstractPartition}

Return whether it is worth to compose `p` and `q`, i.e. the result is less or
equal `max_length` and not equal to the empty partition.

Note that this is a helper function for the `construct_category` algorithm.
"""
function _is_worth_composition(p::T, q::T, max_length::Int) where {T <: AbstractPartition}
    return length(upper_points(p)) != max_length && 
        length(lower_points(p)) + length(upper_points(q)) <= max_length
end

"""
    _add_partition(dict::Dict{Int, Set{T}}, p::T) where {T <: AbstractPartition}

Return `dict` with `p` included as value and `size(p)` as corresponding key.

Note that this is a helper function of the `construct_category` algorithm.
"""
function _add_partition(dict::Dict{Int, Set{T}}, p::T) where 
                                                {T <: AbstractPartition}

    add_apbs = dict[size(p)]
    push!(add_apbs, p)
    dict[size(p)] = add_apbs

    return dict
end

"""
    _add_partition_top_bottom(vector::Vector{Dict{Int, Set{T}}}, p::T) where 
                                                            {T <: AbstractPartition}

Return `vector` with `p` included as value and `size(lower_points(p))` as
corresponding key in the dict of `vector[2]` as well as
`p` included as value and `size(upper_points(p))` as
corresponding key in the dict of `vector[1]`.

Note that this is a helper function of the `construct_category` algorithm.
"""
function _add_partition_top_bottom(vector::Vector{Dict{Int, Set{T}}}, p::T) where 
                                                                {T <: AbstractPartition}

    # add right partition in first dict for top size
    add_apbs_top = (vector[1])[length(upper_points(p))]
    push!(add_apbs_top, p)
    (vector[1])[length(upper_points(p))] = add_apbs_top

    # add right partition in first dict for bottom size
    add_apbs_bottom = (vector[2])[length(upper_points(p))]
    push!(add_apbs_bottom, p)
    (vector[2])[length(lower_points(p))] = add_apbs_bottom

    return vector
end
