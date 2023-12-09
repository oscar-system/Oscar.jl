"""
    helper_new_point_values_vector(p::Vector{Vector{Int}}, q::Vector{Vector{Int}})

Return a semantically identical partition in form of an vector, 
which has new number values.
"""
function helper_new_point_values_vector(p::Vector, q::Vector)

    p_points = vcat(p[1], p[2])

    if !(isempty(p_points))
        new_id, index = findmax(p_points)
        new_id += 1
    else
        new_id = 1 
    end
    
    q_points = vcat(q[1], q[2])
    new_ids = Dict()

    for n in q_points
        new_ids[n] = new_id
        new_id += 1
    end
    upper = [get(new_ids, n, -1) for n in q[1]]
    
    lower = [get(new_ids, n, -1) for n in q[2]]
    
    [upper, lower]
end

"""
    normal_form_vector(p::Vector)

Return a semantically identical partition of vector form, 
which has new number values from 1 to number of blocks of `p`.
"""
function normal_form_vector(p::Vector)

    new_id = 1
    new_ids = Dict{Int, Int}()
    p_return = helper_new_point_values_vector([[length(p[1]) + length(p[2])], []], p)

    for (i, n) in enumerate(p_return[1])
        if !(n in keys(new_ids))
            new_ids[n] = new_id
            p_return[1][i] = new_id
            new_id += 1
        else
            p_return[1][i] = get(new_ids, n, -1)
        end
    end

    for (i, n) in enumerate(p_return[2])
        if !(n in keys(new_ids))
            new_ids[n] = new_id
            p_return[2][i] = new_id
            new_id += 1
        else
            p_return[2][i] = get(new_ids, n, -1)
        end
    end
    p_return
end

"""
    add_partition_to_dict(dict::Dict{Int, Set{AbstractPartition}}, 
        p::AbstractPartition)

Return `dict` with `p` included as value and `size(p)` as
corresponding key.

Note that this is a helper function of the `construct_category` algorithm.
"""
function add_partition_to_dict(dict::Dict{Int, Set{AbstractPartition}}, 
        p::AbstractPartition)

    add_apbs = get(dict, size(p), -1)
    push!(add_apbs, p)
    dict[size(p)] = add_apbs

    dict
end

"""
    add_partition_to_composition_dict(vector::Vector, p::AbstractPartition)

Return `vector` with `p` included as value and `size(lower_points(p))` as
corresponding key in the dict of `vector[2]` as well as
`p` included as value and `size(upper_points(p))` as
corresponding key in the dict of `vector[1]`.

Note that this is a helper function of the `construct_category` algorithm.
"""
function add_partition_to_composition_dict(vector::Vector, p::AbstractPartition)

    # add right partition in first dict for top size
    add_apbs_top = get(vector[1], length(upper_points(p)), -1)
    push!(add_apbs_top, p)
    (vector[1])[length(upper_points(p))] = add_apbs_top

    # add right partition in first dict for bottom size
    add_apbs_bottom = get(vector[2], length(lower_points(p)), -1)
    push!(add_apbs_bottom, p)
    (vector[2])[length(lower_points(p))] = add_apbs_bottom

    vector
end
