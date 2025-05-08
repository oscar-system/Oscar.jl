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
    _is_worth_composition(p::T, q::T, max_length::Int) where {T <: AbstractPartition}

Return whether it is worth to compose `p` and `q`, i.e. the result is less or
equal `max_length` and not equal to the empty partition.

Note that this is a helper function for the `construct_category` algorithm.
"""
function _is_worth_composition(p::T, q::T, max_length::Int) where {T <: AbstractPartition}
    return number_of_upper_points(p) != max_length && 
        number_of_lower_points(p) + number_of_upper_points(q) <= max_length
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
    add_apbs_top = vector[1][number_of_upper_points(p)]
    push!(add_apbs_top, p)
    vector[1][number_of_upper_points(p)] = add_apbs_top

    # add right partition in first dict for bottom size
    add_apbs_bottom = vector[2][number_of_upper_points(p)]
    push!(add_apbs_bottom, p)
    vector[2][number_of_lower_points(p)] = add_apbs_bottom

    return vector
end

"""
    simplify_operation(partition_sum::Vector{Tuple{S, T}}) where { S <: AbstractPartition, T <: RingElement }

Simplify the vector representation of a `LinearPartition` in terms of distributivity.

# Examples
```jldoctest
julia> S, d = polynomial_ring(QQ, :d)
(Univariate polynomial ring in d over QQ, d)

julia> Oscar.SetPartitions.simplify_operation([(set_partition([1, 1], [1, 1]), S(10)), (set_partition([1, 1], [1, 1]), 4*d)])
1-element Vector{Tuple{SetPartition, QQPolyRingElem}}:
 (SetPartition([1, 1], [1, 1]), 4*d + 10)
```
"""
function simplify_operation(partition_sum::Vector{Tuple{S, T}}) where { S <: AbstractPartition, T <: RingElement }

    partitions = Dict{S, T}()

    for (i1, i2) in partition_sum
        if iszero(i2)
            continue
        end
        partitions[i1] = get(partitions, i1, 0) + i2
    end
    
    return [(s, t) for (s, t) in partitions if !iszero(t)]
end

"""
    simplify_operation_zero(p::Dict{S, T}) where { S <: AbstractPartition, T <: RingElement }

Simplify the dict representation of a `LinearPartition` in terms of zero coefficients.
"""
function simplify_operation_zero(p::Dict{S, T}) where { S <: AbstractPartition, T <: RingElement }
    result = Dict{S, T}()
    for (i1, i2) in pairs(p)
        if !iszero(i2)
            result[i1] = i2
        end
    end
    return result
end
