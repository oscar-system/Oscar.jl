# HACK HACK HACK special case optimization...
function GAP.GapObj(obj::Vector{Int}; recursive::Bool = false)
    len = length(obj)
    ret_val = GAP.NewPlist(len)
    for i = 1:len
        ret_val[i] = obj[i]
    end
    return ret_val
end

function normal_form(o::Origami)
    x = horizontal_perm(o)
    y = vertical_perm(o)
    n = degree(o)

    # TODO does this exists already? if not move to appropriate part of Oscar
    function cycle_length(sigma::PermGroupElem, i::Integer)
        return GAPWrap.CYCLE_LENGTH_PERM_INT(GapObj(sigma), i)
    end

    # Find points which minimize the lengths of the cycles in which they occur.
    # This can greatly reduce the number of breadths-first searches below.
    minimal_cycle_lengths = (n, n)
    minimize_cycle_lengths = Int[]
    for i in 1:n
        cycle_lengths = (cycle_length(x, i), cycle_length(y, i))
        if cycle_lengths == minimal_cycle_lengths
            push!(minimize_cycle_lengths, i)
        elseif cycle_lengths < minimal_cycle_lengths
            minimize_cycle_lengths = [i]
            minimal_cycle_lengths = cycle_lengths
        end
    end

    G = PermGroupElem[]
    sym = symmetric_group(n)

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
        push!(G, perm(sym, L))
    end
    
    G2 = [ (x ^ i, y ^ i) for i in G ]
    # no need to check if surface connected
    min_entry = minimum(G2)
    return origami_disconnected(min_entry[1], min_entry[2], n)

end
