function normal_form(o::Origami)
    x = horizontal_perm(o)
    y = vertical_perm(o)
    n = degree(o)

    # TODO does this exists already? if not move to appropriate part of Oscar
    function cycle_length(sigma::PermGroupElem, i::Integer)
        return GAP.Globals.CycleLength(GapObj(sigma), i)::Int
    end

    # Find points which minimize the lengths of the cycles in which they occur.
    # This can greatly reduce the number of breadths-first searches below.
    minimal_cycle_lengths = (n, n)
    minimize_cycle_lengths = Int[]
    degree_list = collect(1:n)
    for i in degree_list
        cycle_lengths = (cycle_length(x, i), cycle_length(y, i))
        if cycle_lengths == minimal_cycle_lengths
            push!(minimize_cycle_lengths, i)
        elseif cycle_lengths < minimal_cycle_lengths
            minimize_cycle_lengths = [i]
            minimal_cycle_lengths = cycle_lengths
        end
    end

    G = []

    # TODO this can be completetly refactored. this is the exact same as
    # normalform_conjugators apart from looping over different lists
    for i in minimize_cycle_lengths
        L = fill(0, n)
        seen = fill(false, n)
        Q = [i]
        seen[i] = true
        numSeen = 1
        L[i] = 1
        while numSeen < n
            v = popfirst!(Q)
            wx = v^x
            wy = v^y
            if !seen[wx]
                push!(Q, wx)
                seen[wx] = true
                numSeen = numSeen + 1
                L[wx] = numSeen
            end
            if !seen[wy]
                push!(Q, wy)
                seen[wy] = true
                numSeen = numSeen + 1
                L[wy] = numSeen
            end
        end
        push!(G, L)
    end
    
    G = perm.(G)
    G = (i -> [i^-1 * x * i, i^-1 * y * i]).(G)
    # no need to check if surface connected
    min_entry = minimum(G)
    return origami_disconnected(min_entry[1], min_entry[2], n)
end