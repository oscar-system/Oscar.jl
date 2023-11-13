"""
    helper_new_point_values_vector(p::vector)

This function outputs a semantically identical Partition in form of an vector which has new number values.

# Arguments
- `p`: The input partition as vector which we do not change
- `q`: The input partition as vector that we change according to the values of p

# Returns
- Semantically equal partition as vector to q without point numbers in p
"""
function helper_new_point_values_vector(p::Vector, q::Vector)

    p_points = vcat(copy(p[1]), copy(p[2]))

    if !(isempty(p_points))
        new_id, index = findmax(p_points)
        new_id += 1
    else
        new_id = 1 
    end
    
    q_points = vcat(copy(q[1]), copy(q[2]))
    new_ids = Dict()

    for n in q_points
        new_ids[n] = new_id
        new_id += 1
    end
    upper = []
    for (i, n) in enumerate(q[1])
        push!(upper, get(new_ids, n, -1))
    end
    
    lower = []
    for (i, n) in enumerate(q[2])
        push!(lower, get(new_ids, n, -1))
    end
    
    [upper, lower]
    
end

"""
    normal_form_vector(p::Vector)

This function outputs a semantically identical Partition of vector form which has new number values from 1 to number of blocks of Partitions.

# Arguments
- `p`: Input partition as vector

# Returns
- `p` with consisten form
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
            p_return[1][i] = get(new_ids, p_return[1][i], -1)
        end
    end

    for (i, n) in enumerate(p_return[2])
        if !(n in keys(new_ids))
            new_ids[n] = new_id
            p_return[2][i] = new_id
            new_id += 1
        else
            p_return[2][i] = get(new_ids, p_return[2][i], -1)
        end
    end
    p_return
end

function add_partition_to_dict(dict::Dict, p::AbstractPartition)

    add_apbs::Set = get(dict, size(p), -1)
    push!(add_apbs, p)
    dict[size(p)] = add_apbs

    return dict
end

function add_partition_to_composition_dict(vector::Vector, p::AbstractPartition)

    # add right partition in first dict for top size
    add_apbs_top::Set = get(vector[1], length(upper_points(p)), -1)
    push!(add_apbs_top, p)
    (vector[1])[length(upper_points(p))] = add_apbs_top

    # add right partition in first dict for bottom size
    add_apbs_bottom::Set = get(vector[2], length(lower_points(p)), -1)
    push!(add_apbs_bottom, p)
    (vector[2])[length(lower_points(p))] = add_apbs_bottom

    return vector
end
