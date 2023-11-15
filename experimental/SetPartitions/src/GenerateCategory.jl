function do_unary(
    to_unary::Set{T}, 
    all_partitions::Set{T}, 
    stop_whole::Bool, 
    already_u::Set{T}, 
    max_length::Int, 
    all_partitions_by_size::Dict{Int, Set{T}}, 
    all_partitions_by_size_top_bottom::Vector{Dict{Int, Set{T}}}, 
    trace::Dict{T, Tuple},
    spatial_rotation::Union{Function, Nothing}) where {T <: AbstractPartition}

    stop = false
    while !stop
        stop = true

        to_unary_copy = copy(to_unary)

        for pp in to_unary_copy
            
            pmod = copy(pp)

            a = copy(pp)

            if pmod isa SetPartition || pmod isa ColoredPartition
                # start with rotation
                if !isempty(upper_points(pmod))
                    a = rotation(pmod, true, true)
                elseif length(lower_points(pmod)) > 0
                    a = rotation(pmod, false, false)
                end

                # add every new partition to all_partitions
                if !(a in all_partitions)
                    trace[a] = tuple([pmod, "r"])
                    stop_whole = false
                    stop = false
                    push!(all_partitions, a)
                    push!(to_unary, a)

                    # call functions which adds partition a into the right set in the dict
                    all_partitions_by_size = 
                        add_partition_to_dict(all_partitions_by_size, a)
                    all_partitions_by_size_top_bottom = add_partition_to_composition_dict(
                            all_partitions_by_size_top_bottom, a)
                end
            elseif spatial_rotation !== nothing
                a = spatial_rotation(pmod)

                # add every new partition to all_partitions
                if !(a in all_partitions)
                    trace[a] = tuple([pmod, "r"])
                    stop_whole = false
                    stop = false
                    push!(all_partitions, a)
                    push!(to_unary, a)

                    # call functions which adds partition a into the right set in the dict
                    all_partitions_by_size = 
                        add_partition_to_dict(all_partitions_by_size, a)
                    all_partitions_by_size_top_bottom = add_partition_to_composition_dict(
                            all_partitions_by_size_top_bottom, a)
                end
            end

            # continue with involution
            a = involution(pmod)

            # add every new partition to all_partitions
            if !(a in all_partitions)
                trace[a] = tuple([pmod, "i"])
                stop_whole = false
                stop = false
                push!(all_partitions, a)
                push!(to_unary, a)

                # call functions which adds the partition a into the right set in the dict
                all_partitions_by_size = add_partition_to_dict(all_partitions_by_size, a)
                all_partitions_by_size_top_bottom = add_partition_to_composition_dict(
                        all_partitions_by_size_top_bottom, a)
            end

            if pmod isa SetPartition
                # end with vertical reflection
                a = vertical_reflection(pmod)

                # add every new partition to all_partitions
                if !(a in all_partitions)
                    trace[a] = tuple([pmod, "vr"])
                    stop_whole = false
                    stop = false
                    push!(all_partitions, a)
                    push!(to_unary, a)

                    # call functions which adds partition a into the right set in the dict 
                    all_partitions_by_size = 
                        add_partition_to_dict(all_partitions_by_size, a)
                    all_partitions_by_size_top_bottom = add_partition_to_composition_dict(
                            all_partitions_by_size_top_bottom, a)
                end
            end
            # remember already processed partitions 
            push!(already_u, pp)
            pop!(to_unary, pp)
        end
    end

    return (stop_whole, all_partitions, already_u, all_partitions_by_size,
        all_partitions_by_size_top_bottom, trace)

end

function do_tensor_products(
    all_partitions::Set{T}, 
    already_t::Set{Tuple}, 
    to_tens::Set{Tuple{T, T}}, 
    stop_whole::Bool, 
    max_length::Int, 
    all_partitions_by_size::Dict{Int, Set{T}}, 
    all_partitions_by_size_top_bottom::Vector{Dict{Int, Set{T}}}, 
    trace::Dict{T, Tuple}) where {T <: AbstractPartition}
    
    # similar to all_pyrtitions_by_size in build function for new_tens
    new_tens_by_size = Dict{Int, Set{AbstractPartition}}(i => Set() for i in 0:max_length)

    # store all partitions which are newly constructed by tensor product
    new_tens = Set{AbstractPartition}()

    # store for every i the ii's which are already used
    without = Dict{AbstractPartition, Set{AbstractPartition}}()

    # until no more new possibilities
    stop = false
    while !stop
        stop = true

        # if there are new partitions due to tensor product and size constraint, 
        # remove pairs which are already calculated
        if !isempty(new_tens)
            aa = union(new_tens, all_partitions)
            for i in aa
                # get fitting partitions in advance (improve runtime)
                new_tens_temp_tensor = Set{AbstractPartition}()
                for key in keys(new_tens_by_size)
                    if size(i) + Int(key) <= max_length
                        new_tens_temp_tensor = union(new_tens_temp_tensor, 
                            get(new_tens_by_size, key, -1))
                    end
                end
                if i in keys(without)
                    operate_on = setdiff(new_tens_temp_tensor, get(without, i, -1))
                else
                    operate_on = new_tens_temp_tensor
                end
                for ii in operate_on
                    if size(i) + size(ii) <= max_length && !((i, ii) in already_t)
                        push!(to_tens, (i, ii))
                        push!(already_t, (i, ii))
                    end
                end
            end
        end

        # do the tensor products
        al = copy(to_tens)
        for (i, ii) in al
            a = tensor_product(i, ii)
            pop!(to_tens, (i, ii))
            if !(a in all_partitions)
                trace[a] = ((i, ii), "t")
                if size(a) == max_length
                    push!(all_partitions, a)

                    # call function which adds partition a to the right set in the dicts
                    all_partitions_by_size = 
                        add_partition_to_dict(all_partitions_by_size, a)
                    all_partitions_by_size_top_bottom = add_partition_to_composition_dict(
                        all_partitions_by_size_top_bottom, a)

                    stop_whole = false
                else
                    push!(all_partitions, a)

                    # call function which adds partition a to the right set in the dicts
                    all_partitions_by_size = 
                        add_partition_to_dict(all_partitions_by_size, a)
                    all_partitions_by_size_top_bottom = add_partition_to_composition_dict(
                        all_partitions_by_size_top_bottom, a)
                    new_tens_by_size = add_partition_to_dict(new_tens_by_size, a)

                    stop_whole = false
                    push!(new_tens, a)
                    stop = false
                end
            else
                # remove not fitting candidates for further iterations
                if !(i in keys(without))
                    without[i] = Set([ii])
                else
                    push!(get(without, i, -1), ii)
                end
            end
        end
    end
    return (all_partitions, already_t, stop_whole, all_partitions_by_size, 
        all_partitions_by_size_top_bottom, trace)
end

function do_composition(
    all_partitions::Set{T}, 
    already_c::Set{Tuple}, 
    stop_whole::Bool, 
    max_length::Int, 
    to_comp::Set{Tuple{T, T}}, 
    all_partitions_by_size::Dict{Int, Set{T}}, 
    all_partitions_by_size_top_bottom::Vector{Dict{Int, Set{T}}}, 
    trace::Dict{T, Tuple}) where {T <: AbstractPartition}

    # add newfound partitions due to composition
    new_comp = Set{AbstractPartition}()

    # new_comp stored in vector with a dict for top and bottom size 
    # (similar to the technique in the construct_category function)
    new_comp_by_size_top_bottom = [Dict{Int, Set{AbstractPartition}}(), Dict{Int, Set{Int}}()]
    new_comp_by_size_top_bottom[1] = Dict{Int, Set{AbstractPartition}}(i => Set() for i in 0:max_length)
    new_comp_by_size_top_bottom[2] = Dict{Int, Set{AbstractPartition}}(i => Set() for i in 0:max_length)
    
    # store for every i the ii's which are already used
    without = Dict{AbstractPartition, Set{AbstractPartition}}()

    # until no more new possibilities
    stop = false
    while !stop
        stop = true
        # if there are new partitions due to composition, 
        # remove pairs which are already calculated
        if !isempty(new_comp)
            aa = union(new_comp, all_partitions)
            for i in aa
                # get fitting partitions in advance (improve runtime)
                new_comp_temp_comp = Set()
                if length(upper_points(i)) <= max_length
                    new_comp_temp_comp = get(new_comp_by_size_top_bottom[2], 
                        length(upper_points(i)), -1)
                end
                if i in keys(without)
                    operate_on = setdiff(new_comp_temp_comp, get(without, i, -1))
                else
                    operate_on = new_comp_temp_comp
                end
                for ii in operate_on
                    if is_composable(i, ii) && is_worth_composition(i, ii, max_length)
                        push!(to_comp, (i, ii))
                        push!(already_c, (i, ii))
                    end
                    if is_composable(ii, i) && is_worth_composition(ii, i, max_length)
                        push!(to_comp, (ii, i))
                        push!(already_c, (ii, i))
                    end
                end
            end
        end

        # do the compositions
        al = copy(to_comp)
        
        for (i, ii) in al
            a = composition(i, ii)
            if !(a in all_partitions)
                trace[a] = ((i, ii), "c")
                push!(all_partitions, a)

                # call function which adds the partition a to the right set in the dicts
                all_partitions_by_size = add_partition_to_dict(all_partitions_by_size, a)
                all_partitions_by_size_top_bottom = add_partition_to_composition_dict(
                    all_partitions_by_size_top_bottom, a)
                new_comp_by_size_top_bottom = 
                    add_partition_to_composition_dict(new_comp_by_size_top_bottom, a)

                stop_whole = false
                push!(new_comp, a)
                stop = false
            else
                # remove not fitting candidates for further iterations
                pop!(to_comp, (i, ii))
                if !(i in keys(without))
                    without[i] = Set([ii])
                else
                    push!(get(without, i, -1), ii)
                end
            end
        end
    end
    return (all_partitions, already_c, stop_whole, all_partitions_by_size, 
        all_partitions_by_size_top_bottom, trace)
end

"""
    construct_category(p::Vector{AbstractPartition}, n::Int, tracing::Bool = false, 
        max_artifical::Int = 0, spatial_rotation::Union{Function,Nothing}=nothing)

Return a list of all partitions of size `n` constructed from partitions in `p` without using partitions of size greater than max(`n`, maxsize(`p`), `max_artifical`)

# Arguments
- `p`: list of partitions
- `n`: size of partitions in constructing category
- `tracing`: optinal input: activate tracing and get the output (category, trace)
- `max_artifical`: optional input: allow partitions to grow > max_length
- `spatial_rotation`: optional input: function from SpatialPartition to SpatialPartition 
    which performs a valid rotation or similar

# Returns
- list of all partitions of size `n` constructed from partitions in `p`

# Examples
```jldoctest
julia> length(construct_category([SetPartition([1, 2], [2, 1])], 6))
105
```
"""
function construct_category(
    p::Vector, 
    n::Int, 
    tracing::Bool = false, 
    max_artificial::Int = 0, 
    spatial_rotation::Union{Function, Nothing} = nothing)
    
    @assert spatial_rotation === nothing || (p isa Vector{SpatialPartition})

    # store all candidates found
    all_partitions = Set{AbstractPartition}()

    # all candidates stored in dict from size to partition
    all_partitions_by_size = Dict{Int, Set{AbstractPartition}}()

    # all candidates stored in vector with a dict for top and bottom size
    all_partitions_by_size_top_bottom = 
        [Dict{Int, Set{AbstractPartition}}(), Dict{Int, Set{AbstractPartition}}()]

    # store partitions already unary
    already_u = Set{AbstractPartition}()

    # store partitions already tensor product
    already_t = Set{Tuple}()

    # store partitions already composition
    already_c = Set{Tuple}()

    # end output: All partitions found of size n
    all_partitions_of_size_n = []

    # all candidates for unary operations
    to_unary = Set{AbstractPartition}(copy(p))

    # trace for tracing
    trace = Dict{AbstractPartition, Tuple}()

    # compare allowed expansion size with max(n, max_length)
    max_length = max(n, maximum(size(i) for i in p))

    if max_artificial > 0
        max_length = max(max_length, max_artificial)
    end

    # define for all i <= size an empty set,
    # in which we fill the corresponding partition of size i (for tensor)
    for i in 0:max_length
        all_partitions_by_size[i] = Set{AbstractPartition}()
    end

    # for all sizes bottom and top empty set in which we fill the corresponding partition
    all_partitions_by_size_top_bottom[1] = 
        Dict{Int, Set{AbstractPartition}}(i => Set() for i in 0:max_length)
    all_partitions_by_size_top_bottom[2] = 
        Dict{Int, Set{AbstractPartition}}(i => Set() for i in 0:max_length)

    # add partitions in p to all_partitions_by_size and all_partitions_by_size_top_bottom
    tuple_list_all_partitions = []
    for i in all_partitions
        push!(tuple_list_all_partitions, i)
    end
    for i in vcat(p, tuple_list_all_partitions)
        all_partitions_by_size = add_partition_to_dict(all_partitions_by_size, i)
        all_partitions_by_size_top_bottom = 
            add_partition_to_composition_dict(all_partitions_by_size_top_bottom, i)
    end
    
    # while new were found apply on them unary tensor and composition
    stop_whole = false
    while !stop_whole
        
        stop_whole = true
        # add new found partitions in the unary operation candidate list
        for i in all_partitions
            if !(i in already_u)
                push!(to_unary, i)
            end
        end

        # fist phase: all possible combinations of unary operations
        (stop_whole, all_partitions, already_u, all_partitions_by_size, 
            all_partitions_by_size_top_bottom, trace) = 
            do_unary(to_unary, all_partitions, stop_whole, already_u, max_length, 
            all_partitions_by_size, all_partitions_by_size_top_bottom, trace, 
            spatial_rotation)

        # store pairs as candidates for the tensor products
        to_tens = Set{Tuple{AbstractPartition, AbstractPartition}}()

        # get all candidates
        for i in all_partitions
            # get fitting partitions in advance (improve runtime)
            all_partitions_temp_tensor = Set()
            for key in keys(all_partitions_by_size)
                if size(i) + Int(key) <= max_length
                    all_partitions_temp_tensor = union(all_partitions_temp_tensor, 
                        get(all_partitions_by_size, key, -1))
                end
            end
            for ii in all_partitions_temp_tensor
                if (!((i, ii) in already_t) && (length(upper_points(i)) + 
                    length(lower_points(i)) + length(upper_points(ii)) + 
                    length(lower_points(ii)) <= max_length))
                    push!(to_tens, (i, ii))
                    push!(already_t, (i, ii))
                end
            end
        end
        
        # second phase: all possible tensor products which aren't redundant
        (all_partitions, already_t, stop_whole, all_partitions_by_size, 
            all_partitions_by_size_top_bottom, trace) = 
            do_tensor_products(all_partitions, already_t, to_tens, stop_whole, max_length, 
            all_partitions_by_size, all_partitions_by_size_top_bottom, trace)

        # add new variations produced by tensor product for composition
        to_comp = Set{Tuple{AbstractPartition, AbstractPartition}}()

        # get all candidates for
        for i in all_partitions
            # get in advance the right second candidate (regarding format)
            all_partitions_temp_comp = get(all_partitions_by_size_top_bottom[2], 
                length(upper_points(i)), -1)
            for ii in all_partitions_temp_comp
                if !((i, ii) in already_c) && is_composable(i, ii) && 
                        is_worth_composition(i, ii, max_length)
                    push!(to_comp, (i, ii))
                    push!(already_c, (i, ii))
                end
            end
        end

        # third phase: all possible compositions which aren't redundant
        (all_partitions, already_c, stop_whole, all_partitions_by_size, 
            all_partitions_by_size_top_bottom, trace) = do_composition(all_partitions, 
            already_c, stop_whole, max_length, to_comp, all_partitions_by_size, 
            all_partitions_by_size_top_bottom, trace)
    end

    # remove all partitions without size n
    for i in all_partitions
        if size(i) == n && !(i in all_partitions_of_size_n)
            push!(all_partitions_of_size_n, i)
        end
    end

    if tracing
        return all_partitions_of_size_n, trace
    end

    return all_partitions_of_size_n
end

"""
    get_trace(trace::Dict, start::Partition)

Print the trace of the partition `start` constructed with `construct_category`.
(via breath first search)

"""
function get_trace(trace::Dict{AbstractPartition, Tuple}, start)
    # iterate through trace with breath first search and print it

    if !(start in keys(trace))
        print("(spatial) Partition $(start) not found in trace")
    end
    
    track = [start]
    for p in track
        if p in keys(trace)
            println(p, " : ", get(trace, p, -1))
            for i in get(trace, p, -1)[1]
                if !(typeof(i) <: AbstractString) && !(i in track)
                    push!(track, i)
                end
            end
        end
    end
end
