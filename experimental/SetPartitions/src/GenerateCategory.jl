"""
    Trace{T}

Returned by `construct_category` to trace the construction of partitions.

Each partition of type `T` is mapped to to a `Vector` of partitions from which 
it is constructed and a `String` describing the corresponding operation.
"""
const Trace{T} = Dict{T, Tuple{Vector{T}, String}}

"""
This function is a helper function of `construct_category`, which applies unary operations
(i.e. rotations, involutions, vertical reflections) on all partitions in 
the set `to_unary`. When a new partition is created in the course of this process, 
it is incorporated into `all_partitions`, `all_partitions_by_size`, and 
`all_partitions_by_size_top_bottom`, simultaneously setting the `stop_whole` flag to true.

The return value of this function is a tuple of the updated sets and the `stop_whole` flag.
"""
function _do_unary(
    to_unary::Set{T}, 
    all_partitions::Set{T}, 
    stop_whole::Bool, 
    already_u::Set{T}, 
    max_length::Int, 
    all_partitions_by_size::Dict{Int, Set{T}}, 
    all_partitions_by_size_top_bottom::Vector{Dict{Int, Set{T}}}, 
    trace::Trace{T},
    spatial_rotation::Union{Function, Nothing}) where {T <: AbstractPartition}

    stop = false
    while !stop
        stop = true

        for pp in to_unary
            
            pmod = deepcopy(pp)

            a = deepcopy(pp)

            if pmod isa SetPartition || pmod isa ColoredPartition
                # start with rotation
                if !isempty(upper_points(pmod))
                    a = rotate_top_left(pmod)
                elseif length(lower_points(pmod)) > 0
                    a = rotate_bottom_right(pmod)
                end

                # add every new partition to all_partitions
                if !(a in all_partitions)
                    trace[a] = ([pmod], "r")
                    stop_whole = false
                    stop = false
                    push!(all_partitions, a)
                    push!(to_unary, a)

                    # call functions which adds partition a into the right set in the dict
                    all_partitions_by_size = 
                        _add_partition(all_partitions_by_size, a)
                    all_partitions_by_size_top_bottom = _add_partition_top_bottom(
                            all_partitions_by_size_top_bottom, a)
                end
            elseif spatial_rotation !== nothing
                a = spatial_rotation(pmod)

                # add every new partition to all_partitions
                if !(a in all_partitions)
                    trace[a] = ([pmod], "r")
                    stop_whole = false
                    stop = false
                    push!(all_partitions, a)
                    push!(to_unary, a)

                    # call functions which adds partition a into the right set in the dict
                    all_partitions_by_size = 
                        _add_partition(all_partitions_by_size, a)
                    all_partitions_by_size_top_bottom = _add_partition_top_bottom(
                            all_partitions_by_size_top_bottom, a)
                end
            end

            # continue with involution
            a = involution(pmod)

            # add every new partition to all_partitions
            if !(a in all_partitions)
                trace[a] = ([pmod], "i")
                stop_whole = false
                stop = false
                push!(all_partitions, a)
                push!(to_unary, a)

                # call functions which adds the partition a into the right set in the dict
                all_partitions_by_size = _add_partition(all_partitions_by_size, a)
                all_partitions_by_size_top_bottom = _add_partition_top_bottom(
                        all_partitions_by_size_top_bottom, a)
            end

            if pmod isa SetPartition
                # end with vertical reflection
                a = reflect_vertical(pmod)

                # add every new partition to all_partitions
                if !(a in all_partitions)
                    trace[a] = ([pmod], "vr")
                    stop_whole = false
                    stop = false
                    push!(all_partitions, a)
                    push!(to_unary, a)

                    # call functions which adds partition a into the right set in the dict 
                    all_partitions_by_size = 
                        _add_partition(all_partitions_by_size, a)
                    all_partitions_by_size_top_bottom = _add_partition_top_bottom(
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

"""
This function is a helper function of `construct_category`, which applies 
tensor products on all partitions in the set `to_tens`. When a new partition 
is created in the course of this process, it is incorporated into `all_partitions`, 
`all_partitions_by_size`, `already_t` and `all_partitions_by_size_top_bottom`, 
simultaneously setting the `stop_whole` flag to true.

The return value of this function is a tuple of the updated sets and the `stop_whole` flag.
"""
function _do_tensor_products(
    all_partitions::Set{T}, 
    already_t::Set{Tuple{T, T}}, 
    to_tens::Set{Tuple{T, T}}, 
    stop_whole::Bool, 
    max_length::Int, 
    all_partitions_by_size::Dict{Int, Set{T}}, 
    all_partitions_by_size_top_bottom::Vector{Dict{Int, Set{T}}}, 
    trace::Trace{T}) where {T <: AbstractPartition}
    
    # similar to all_pyrtitions_by_size in build function for new_tens
    new_tens_by_size = Dict{Int, Set{T}}(i => Set{T}() for i in 0:max_length)

    # store all partitions which are newly constructed by tensor product
    new_tens = Set{T}()

    # store for every i the ii's which are already used
    without = Dict{T, Set{T}}()

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
                new_tens_temp_tensor = Set{T}()
                for key in keys(new_tens_by_size)
                    if size(i) + Int(key) <= max_length
                        new_tens_temp_tensor = union(new_tens_temp_tensor, 
                            new_tens_by_size[key])
                    end
                end
                if i in keys(without)
                    operate_on = setdiff(new_tens_temp_tensor, without[i])
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
        al = deepcopy(to_tens)
        for (i, ii) in al
            a = tensor_product(i, ii)
            pop!(to_tens, (i, ii))
            if !(a in all_partitions)
                trace[a] = ([i, ii], "t")
                if size(a) == max_length
                    push!(all_partitions, a)

                    # call function which adds partition a to the right set in the dicts
                    all_partitions_by_size = 
                        _add_partition(all_partitions_by_size, a)
                    all_partitions_by_size_top_bottom = _add_partition_top_bottom(
                        all_partitions_by_size_top_bottom, a)

                    stop_whole = false
                else
                    push!(all_partitions, a)

                    # call function which adds partition a to the right set in the dicts
                    all_partitions_by_size = 
                        _add_partition(all_partitions_by_size, a)
                    all_partitions_by_size_top_bottom = _add_partition_top_bottom(
                        all_partitions_by_size_top_bottom, a)
                    new_tens_by_size = _add_partition(new_tens_by_size, a)

                    stop_whole = false
                    push!(new_tens, a)
                    stop = false
                end
            else
                # remove not fitting candidates for further iterations
                if !(i in keys(without))
                    without[i] = Set{T}([ii])
                else
                    push!(without[i], ii)
                end
            end
        end
    end
    return (all_partitions, already_t, stop_whole, all_partitions_by_size, 
        all_partitions_by_size_top_bottom, trace)
end


"""
This function is a helper function of `construct_category`, which applies 
compositions on all partitions in the set `to_comp`. When a new partition 
is created in the course of this process, it is incorporated into `all_partitions`, 
`all_partitions_by_size`, `already_c` and `all_partitions_by_size_top_bottom`, 
simultaneously setting the `stop_whole` flag to true.

The return value of this function is a tuple of the updated sets and the `stop_whole` flag.
"""
function _do_composition(
    all_partitions::Set{T}, 
    already_c::Set{Tuple{T, T}}, 
    stop_whole::Bool, 
    max_length::Int, 
    to_comp::Set{Tuple{T, T}}, 
    all_partitions_by_size::Dict{Int, Set{T}}, 
    all_partitions_by_size_top_bottom::Vector{Dict{Int, Set{T}}}, 
    trace::Trace{T}) where {T <: AbstractPartition}

    # add newfound partitions due to composition
    new_comp = Set{T}()

    # new_comp stored in vector with a dict for top and bottom size 
    # (similar to the technique in the construct_category function)
    new_comp_by_size_top_bottom = 
        [Dict{Int, Set{T}}(), Dict{Int, Set{T}}()]
    new_comp_by_size_top_bottom[1] = 
        Dict{Int, Set{T}}(i => Set{T}() for i in 0:max_length)
    new_comp_by_size_top_bottom[2] = 
        Dict{Int, Set{T}}(i => Set{T}() for i in 0:max_length)
    
    # store for every i the ii's which are already used
    without = Dict{T, Set{T}}()

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
                new_comp_temp_comp = Set{T}()
                if length(upper_points(i)) <= max_length
                    new_comp_temp_comp = 
                        (new_comp_by_size_top_bottom[2])[length(upper_points(i))]
                end
                if i in keys(without)
                    operate_on = setdiff(new_comp_temp_comp, without[i])
                else
                    operate_on = new_comp_temp_comp
                end
                for ii in operate_on
                    if is_composable(i, ii) && _is_worth_composition(i, ii, max_length)
                        push!(to_comp, (i, ii))
                        push!(already_c, (i, ii))
                    end
                    if is_composable(ii, i) && _is_worth_composition(ii, i, max_length)
                        push!(to_comp, (ii, i))
                        push!(already_c, (ii, i))
                    end
                end
            end
        end

        # do the compositions
        al = deepcopy(to_comp)
        
        for (i, ii) in al
            a = compose(i, ii)
            if !(a in all_partitions)
                trace[a] = ([i, ii], "c")
                push!(all_partitions, a)

                # call function which adds the partition a to the right set in the dicts
                all_partitions_by_size = _add_partition(all_partitions_by_size, a)
                all_partitions_by_size_top_bottom = _add_partition_top_bottom(
                    all_partitions_by_size_top_bottom, a)
                new_comp_by_size_top_bottom = 
                    _add_partition_top_bottom(new_comp_by_size_top_bottom, a)

                stop_whole = false
                push!(new_comp, a)
                stop = false
            else
                # remove not fitting candidates for further iterations
                pop!(to_comp, (i, ii))
                if !(i in keys(without))
                    without[i] = Set{T}([ii])
                else
                    push!(without[i], ii)
                end
            end
        end
    end
    return (all_partitions, already_c, stop_whole, all_partitions_by_size, 
        all_partitions_by_size_top_bottom, trace)
end

"""
    construct_category(p::Vector{AbstractPartition}, n::Int, tracing::Bool = false, 
        max_artifical::Int = 0, spatial_rotation::Union{Functi on,Nothing}=nothing)

Return a list of all partitions of size `n` which can be constructed from category 
operations using partitions in `p` and without using partitions of size greater than 
`max(n, maxsize(p), max_artifical)`.

Category operations include composition, tensor product, involution, 
rotation and reflection. See Section 4.1.1 in [Gro20](@cite) 
for more information on categories of partitions and these operations. 

See also Section 4 in [Vol23](@cite) for a description of the underlying algorithm.

# Arguments
- `p`: list of partitions
- `n`: size of partitions to construct
- `tracing` (optional): return additional data to allow tracing using `print_trace` 
- `max_artifical` (optional): allow partitions to grow larger then `n` and `maxsize(p)`
- `spatial_rotation` (optional): function which performs a rotation on `SpatialPartition`

# Returns
- list of all partitions of size `n` constructed from partitions in `p`

# Examples
```jldoctest
julia> length(construct_category([SetPartition([1, 2], [2, 1])], 6))
105
```
"""
function construct_category(
    p::Vector{T}, 
    n::Int, 
    tracing::Bool = false, 
    max_artificial::Int = 0, 
    spatial_rotation::Union{Function, Nothing} = nothing) where {T <: AbstractPartition}
    
    @assert spatial_rotation === nothing || (p isa Vector{SpatialPartition})

    # store all candidates found
    all_partitions = Set{T}()

    # all candidates stored in dict from size to partition
    all_partitions_by_size = Dict{Int, Set{T}}()

    # all candidates stored in vector with a dict for top and bottom size
    all_partitions_by_size_top_bottom = 
        [Dict{Int, Set{T}}(), Dict{Int, Set{T}}()]

    # store partitions already unary
    already_u = Set{T}()

    # store partitions already tensor product
    already_t = Set{Tuple{T, T}}()

    # store partitions already composition
    already_c = Set{Tuple{T, T}}()

    # end output: All partitions found of size n
    all_partitions_of_size_n = []

    # all candidates for unary operations
    to_unary = Set{T}(deepcopy(p))

    # trace for tracing
    trace = Trace{T}()

    # compare allowed expansion size with max(n, max_length)
    max_length = max(n, maximum(size(i) for i in p))

    if max_artificial > 0
        max_length = max(max_length, max_artificial)
    end

    # define for all i <= size an empty set,
    # in which we fill the corresponding partition of size i (for tensor)
    for i in 0:max_length
        all_partitions_by_size[i] = Set{T}()
    end

    # for all sizes bottom and top empty set in which we fill the corresponding partition
    all_partitions_by_size_top_bottom[1] = 
        Dict{Int, Set{T}}(i => Set{T}() for i in 0:max_length)
    all_partitions_by_size_top_bottom[2] = 
        Dict{Int, Set{T}}(i => Set{T}() for i in 0:max_length)

    # add partitions in p to all_partitions_by_size and all_partitions_by_size_top_bottom
    tuple_list_all_partitions = []
    for i in all_partitions
        push!(tuple_list_all_partitions, i)
    end
    for i in vcat(p, tuple_list_all_partitions)
        all_partitions_by_size = _add_partition(all_partitions_by_size, i)
        all_partitions_by_size_top_bottom = 
            _add_partition_top_bottom(all_partitions_by_size_top_bottom, i)
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
            _do_unary(to_unary, all_partitions, stop_whole, already_u, max_length, 
            all_partitions_by_size, all_partitions_by_size_top_bottom, trace, 
            spatial_rotation)

        # store pairs as candidates for the tensor products
        to_tens = Set{Tuple{T, T}}()

        # get all candidates
        for i in all_partitions
            # get fitting partitions in advance (improve runtime)
            all_partitions_temp_tensor = Set{T}()
            for key in keys(all_partitions_by_size)
                if size(i) + Int(key) <= max_length
                    all_partitions_temp_tensor = union(all_partitions_temp_tensor, 
                        all_partitions_by_size[key])
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
            _do_tensor_products(all_partitions, already_t, to_tens, stop_whole, max_length,
            all_partitions_by_size, all_partitions_by_size_top_bottom, trace)

        # add new variations produced by tensor product for composition
        to_comp = Set{Tuple{T, T}}()

        # get all candidates for
        for i in all_partitions
            # get in advance the right second candidate (regarding format)
            all_partitions_temp_comp = 
                (all_partitions_by_size_top_bottom[2])[length(upper_points(i))]
            for ii in all_partitions_temp_comp
                if !((i, ii) in already_c) && is_composable(i, ii) && 
                        _is_worth_composition(i, ii, max_length)
                    push!(to_comp, (i, ii))
                    push!(already_c, (i, ii))
                end
            end
        end

        # third phase: all possible compositions which aren't redundant
        (all_partitions, already_c, stop_whole, all_partitions_by_size, 
            all_partitions_by_size_top_bottom, trace) = _do_composition(all_partitions, 
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
    print_trace(trace::Trace{T}, start::T) where {T <: AbstractPartition}

Print the trace of the partition `start` constructed with `construct_category`.
"""
function print_trace(trace::Trace{T}, start::T) where {T <: AbstractPartition}
    # iterate through trace with breath first search and print it

    if !(start in keys(trace))
        print("(spatial) Partition $(start) not found in trace")
    end
    
    track = [start]
    for p in track
        if p in keys(trace)
            println(p, " : ", trace[p])
            for i in (trance[p])[1]
                if !(i in track)
                    push!(track, i)
                end
            end
        end
    end
end
