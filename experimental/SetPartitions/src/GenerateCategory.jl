function do_unary(
    to_unary::Set, 
    all_partitions::Set, 
    stop_whole::Bool, 
    already_u::Set, 
    max_length::Int, 
    all_partitions_by_size::Dict, 
    all_partitions_by_size_top_bottom::Array, 
    trace::Dict,
    spatial_rotation::Function)

    stop::Bool = false
    while !stop
        stop = true

        to_unary_copy::Set{AbstractPartition} = copy(to_unary)

        for pp in to_unary_copy
            
            pmod::AbstractPartition = copy(pp)

            a::AbstractPartition = copy(pp)

            if pmod isa Partition || pmod isa ColoredPartition
                # start with rotation
                if !isempty(get_upper_points(pmod))
                    a = rotation(pmod, true, true)
                elseif length(get_lower_points(pmod)) > 0
                    a = rotation(pmod, false, false)
                end

                # add to all_partitions
                if !(a in all_partitions)
                    trace[a] = tuple([pmod, "r"])
                    stop_whole = false
                    stop = false
                    push!(all_partitions, a)
                    push!(to_unary, a)

                    # call functions which adds the partition a into the right set in the dict
                    all_partitions_by_size = add_partition_to_dict(all_partitions_by_size, a)
                    all_partitions_by_size_top_bottom = add_partition_to_composition_dict(
                            all_partitions_by_size_top_bottom, a)
                end
            elseif length(methods(spatial_rotation)[1].sig.parameters)-1 == 1
                a = spatial_rotation(pmod)

                # add to all_partitions
                if !(a in all_partitions)
                    trace[a] = tuple([pmod, "r"])
                    stop_whole = false
                    stop = false
                    push!(all_partitions, a)
                    push!(to_unary, a)

                    # call functions which adds the partition a into the right set in the dict
                    all_partitions_by_size = add_partition_to_dict(all_partitions_by_size, a)
                    all_partitions_by_size_top_bottom = add_partition_to_composition_dict(
                            all_partitions_by_size_top_bottom, a)
                end
            end

            # continue with involution
            a = involution(pmod)

            # add to all_partitions
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

            if pmod isa Partition
                # end with involution y-axis
                a = vertical_reflection(pmod)

                # add to all_partitions
                if !(a in all_partitions)
                    trace[a] = tuple([pmod, "vr"])
                    stop_whole = false
                    stop = false
                    push!(all_partitions, a)
                    push!(to_unary, a)

                    # call functions which adds the partition a into the right set in the dict 
                    all_partitions_by_size = add_partition_to_dict(all_partitions_by_size, a)
                    all_partitions_by_size_top_bottom = add_partition_to_composition_dict(
                            all_partitions_by_size_top_bottom, a)
                end
            end
            # remember already unary
            push!(already_u, pp)
            pop!(to_unary, pp)
        end
    end

    return [stop_whole, all_partitions, already_u, all_partitions_by_size, all_partitions_by_size_top_bottom, trace]

end

function do_tensor_products(
    all_partitions::Set, 
    already_t::Set, 
    to_tens::Set, 
    stop_whole::Bool, 
    max_length::Int, 
    all_partitions_by_size::Dict, 
    all_partitions_by_size_top_bottom::Array, 
    trace::Dict)::Array
    
    # analogical to all_pyrtitions_by_size in build function for new_tens
    new_tens_by_size = Dict()
    new_tens_by_size = Dict(i => Set() for i in 0:max_length)

    # store all partitions which are new constructed by tensor product
    new_tens::Set{AbstractPartition} = Set()

    # store for every i the ii's which are already used, to not use them in this iteration again
    without::Dict{AbstractPartition, Set{AbstractPartition}} = Dict()

    # until no more new possibilities tensor
    stop::Bool = false
    while !stop
        stop = true

        # if there are new partitions due to tensor and size constraint, remove pair which are already 
        # calculated
        if !isempty(new_tens)
            aa::Set{AbstractPartition} = union(new_tens, all_partitions)
            for i in aa
                # get fitting partitions in advance (improve runtime)
                new_tens_temp_tensor::Set{AbstractPartition} = Set()
                for key in keys(new_tens_by_size)
                    if size(i) + Int(key) <= max_length
                        new_tens_temp_tensor = union(new_tens_temp_tensor, get(new_tens_by_size, key, -1)::Set)
                    end
                end
                if i in keys(without)
                    operate_on = setdiff(new_tens_temp_tensor, get(without, i, -1))
                else
                    operate_on = new_tens_temp_tensor
                end
                for ii in operate_on
                    if size(i) + size(ii) <= max_length && !([i, ii] in already_t)
                        push!(to_tens, [i, ii])
                        push!(already_t, [i, ii])
                    end
                end
            end
        end

        # do the tensor products
        al::Set = copy(to_tens)
        for (i, ii) in al
            a::AbstractPartition = tensor_product(i, ii)
            pop!(to_tens, [i, ii])
            if !(a in all_partitions)
                trace[a] = ((i, ii), "t")
                if size(a) == max_length
                    push!(all_partitions, a)

                    # call function which adds the partition a into the right set in the dicts
                    all_partitions_by_size = add_partition_to_dict(all_partitions_by_size, a)
                    all_partitions_by_size_top_bottom = add_partition_to_composition_dict(
                        all_partitions_by_size_top_bottom, a)

                    stop_whole = false
                else
                    push!(all_partitions, a)

                    # call function which adds the partition a into the right set in the dicts
                    all_partitions_by_size = add_partition_to_dict(all_partitions_by_size, a)
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
    return [all_partitions, already_t, stop_whole, all_partitions_by_size, all_partitions_by_size_top_bottom, trace]
end

function do_composition(
    all_partitions::Set, 
    already_c::Set, 
    stop_whole::Bool, 
    max_length::Int, 
    to_comp::Set, 
    all_partitions_by_size::Dict, 
    all_partitions_by_size_top_bottom::Array, 
    trace::Dict)::Array

    # add newfound partitions due comp
    new_comp::Set{AbstractPartition} = Set()

    # new_comp stored in tuple with a dict for top and bottom size (analogical to the technique in build function)
    new_comp_by_size_top_bottom = [Dict(), Dict()]
    new_comp_by_size_top_bottom[1] = Dict(i => Set() for i in 0:max_length)
    new_comp_by_size_top_bottom[2] = Dict(i => Set() for i in 0:max_length)
    
    # store for every i the ii's which are already used, to not use them in this iteration again
    without::Dict{AbstractPartition, Set{AbstractPartition}} = Dict()

    # until no more new possibilities compose
    stop::Bool = false
    while !stop
        #println("while comp", new_comp)
        stop = true
        # if there are new partitions due to composition, remove pair which are already calculated
        if !isempty(new_comp)
            aa = union(new_comp, all_partitions)
            for i in aa
                # get fitting partitions in advance (improve runtime)
                new_comp_temp_comp = Set()
                if length(get_upper_points(i)) <= max_length
                    new_comp_temp_comp = get(new_comp_by_size_top_bottom[2], length(get_upper_points(i)), -1)
                end
                if i in keys(without)
                    operate_on = setdiff(new_comp_temp_comp, get(without, i, -1))
                else
                    operate_on = new_comp_temp_comp
                end
                for ii in operate_on
                    if is_composable(i, ii) && is_worth_composition(i, ii, max_length)
                        push!(to_comp, [i, ii])
                        push!(already_c, [i, ii])
                    end
                    if is_composable(ii, i) && is_worth_composition(ii, i, max_length)
                        push!(to_comp, [ii, i])
                        push!(already_c, [ii, i])
                    end
                end
            end
        end

        # do the compositions
        al = copy(to_comp)

        #println(to_comp)
        
        for (i, ii) in al
            a = composition(i, ii)
            if !(a in all_partitions)
                trace[a] = ((i, ii), "c")
                push!(all_partitions, a)

                # call function which adds the partition a into the right set in the dicts
                all_partitions_by_size = add_partition_to_dict(all_partitions_by_size, a)
                all_partitions_by_size_top_bottom = add_partition_to_composition_dict(
                    all_partitions_by_size_top_bottom, a)
                new_comp_by_size_top_bottom = add_partition_to_composition_dict(new_comp_by_size_top_bottom, a)

                stop_whole = false
                push!(new_comp, a)
                stop = false
            else
                # remove not fitting candidates for further iterations
                pop!(to_comp, [i, ii])
                if !(i in keys(without))
                    without[i] = Set([ii])
                else
                    push!(get(without, i, -1), ii)
                end
            end
        end
    end
    return [all_partitions, already_c, stop_whole, all_partitions_by_size, all_partitions_by_size_top_bottom, trace]
end

"""
construct_category(p::Array, n::Int, tracing::Bool = false, max_artifical::Int = 0)

This function outputs list of all partitions size n constructed from partitions in p (without using partitions of size
greater than max(max(n, max(p)), max_artifical))

# Arguments
- `p`: list of partitions
- `n`: size of partitions in constructing category
- `tracing`: optinal input: activate tracing and get the output (category, trace)
- `max_artifical`: optional input: allow partitions to grow greater while construction process

# Returns
- list of all partitions size n constructed from partitions in p

# Examples
```julia-repl
julia> length(construct_category([Partition([1, 2], [2, 1])], 6))
105
julia> length(construct_category([Partition([1, 2], [2, 1])], 6, true))
[<Partition([1, 2], [2, 1])> ∩ P(6), Dict{Partition -> Tuple}]
```
"""
function construct_category(p::Array, n::Int, tracing::Bool = false, max_artificial::Int = 0, spatial_rotation::Function = nothing_function() = nothing)
    
    @assert length(methods(spatial_rotation)[1].sig.parameters)-1 == 0 || (p isa Array{SpatialPartition})

    # store all candidates found
    all_partitions::Set{AbstractPartition} = Set()

    # all candidates stored in dict from size to partition
    all_partitions_by_size::Dict{Int, Set{AbstractPartition}} = Dict()

    # all candidates stored in tuple with a dict for top and bottom size
    all_partitions_by_size_top_bottom::Array = [Dict(), Dict()]

    # store partitions already unary
    already_u::Set{AbstractPartition} = Set()

    # store partitions already tensor product
    already_t::Set{Array{AbstractPartition}} = Set()

    # store partitions already composition
    already_c::Set{Array{AbstractPartition}} = Set()

    # end output: All partitions found of size n
    all_partitions_of_size_n::Array{AbstractPartition} = []

    # all candidates for unary operations
    to_unary::Set{AbstractPartition} = Set(copy(p))

    # trace for tracing
    trace::Dict = Dict()

    # compare allowed expansion size with max(n, max_length)
    max_length::Int = max(n, maximum(size(i) for i in p))

    if max_artificial > 0
        max_length = max(max_length, max_artificial)
    end

    # define for all i <= size an empty set in which we fill the corresponding partition of size i (for tensor)
    for i in 0:max_length
        all_partitions_by_size[i] = Set()
    end

    # define for all bottom and top size an empty set in which we fill the corresponding partition
    all_partitions_by_size_top_bottom[1] = Dict(i => Set() for i in 0:max_length)
    all_partitions_by_size_top_bottom[2] = Dict(i => Set() for i in 0:max_length)

    # add all partitions in p to all_partitions_by_size and all_partitions_by_size_top_bottom
    tuple_list_all_partitions::Array = []
    for i in all_partitions
        push!(tuple_list_all_partitions, i)
    end
    for i in vcat(p, tuple_list_all_partitions)
        all_partitions_by_size = add_partition_to_dict(all_partitions_by_size, i)
        all_partitions_by_size_top_bottom = add_partition_to_composition_dict(all_partitions_by_size_top_bottom, i)
    end
    
    # while new were found apply on them unary tensor and composition
    stop_whole::Bool = false
    while !stop_whole
        
        stop_whole = true
        # add new found partitions in the unary operation candidate list
        for i in all_partitions
            if !(i in already_u)
                push!(to_unary, i)
            end
        end

        # fist phase: all possible combinations of unary operations
        (stop_whole, all_partitions, already_u, all_partitions_by_size, all_partitions_by_size_top_bottom, trace) = do_unary(to_unary, all_partitions, stop_whole, already_u, max_length, all_partitions_by_size, all_partitions_by_size_top_bottom, trace, spatial_rotation)

        # store pairs that are candidates to get tensor product
        to_tens::Set{Array} = Set()

        # get all pairs to tensor
        for i in all_partitions
            # get fitting partitions in advance (improve runtime)
            all_partitions_temp_tensor = Set()
            for key in keys(all_partitions_by_size)
                if size(i) + Int(key) <= max_length
                    all_partitions_temp_tensor = union(all_partitions_temp_tensor, get(all_partitions_by_size, key, -1))
                end
            end
            for ii in all_partitions_temp_tensor
                if !([i, ii] in already_t) && (length(get_upper_points(i)) + length(get_lower_points(i)) + length(get_upper_points(ii)) + length(get_lower_points(ii)) <= max_length)
                    push!(to_tens, [i, ii])
                    push!(already_t, [i, ii])
                end
            end
        end
        
        # second phase: all possible tensor product operations which aren't redundant (don't do tensor products twice)
        (all_partitions, already_t, stop_whole, all_partitions_by_size, all_partitions_by_size_top_bottom, trace) = do_tensor_products(all_partitions, already_t, to_tens, stop_whole, max_length, all_partitions_by_size, all_partitions_by_size_top_bottom, trace)

        # add new variations by tensor product or composition with all others
        to_comp::Set{Array} = Set()

        # get all pairs to compose
        for i in all_partitions
            # get in advance the right second candidate (regarding format)
            all_partitions_temp_comp = get(all_partitions_by_size_top_bottom[2], length(get_upper_points(i)), -1)
            for ii in all_partitions_temp_comp
                if !([i, ii] in already_c) && is_composable(i, ii) && is_worth_composition(i, ii, max_length)
                    push!(to_comp, [i, ii])
                    push!(already_c, [i, ii])
                end
            end
        end

        # third phase: all possible compositions which aren't redundant (don't do tensor products twice)
        (all_partitions, already_c, stop_whole, all_partitions_by_size, all_partitions_by_size_top_bottom, trace) = do_composition(all_partitions, already_c, stop_whole, max_length, to_comp, all_partitions_by_size, all_partitions_by_size_top_bottom, trace)
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

This function prints out the trace of the partition `start` constructed with `construct_category`
(via breath first search)

"""
function get_trace(trace::Dict, start)
    # track the trace with breath first search

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