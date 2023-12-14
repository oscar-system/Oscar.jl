function in_deck_group(o::Origami, sigma::PermGroupElem)
    degree_list = collect(1:degree(o))
    h = horizontal_perm(o)
    v = vertical_perm(o)
    for i in degree_list
        if (i^sigma)^h != (i^h)^sigma || (i^sigma)^v != (i^v)^sigma
            return false
        end
    end
    return true
end

function deck_group(o::Origami)
    function candidate_for_deck(origami, j)
        sheets_to_visit = collect(2:degree(origami))
        sigma = [j]
        append!(sigma, fill(-1, degree(origami) - 1))
        while length(sheets_to_visit) > 0
            found_predecessor = false
            h = horizontal_perm(origami)
            v = vertical_perm(origami)
            for i in sheets_to_visit
                for tau in [h, v]
                    position_iterator = searchsorted(sheets_to_visit, i^(tau^-1))
                    if first(position_iterator) > last(position_iterator)
                        found_predecessor = true
                        sigma[i] = sigma[i^(tau^-1)]^tau
                        deleteat!(sheets_to_visit, findfirst(item -> item == i, sheets_to_visit))
                        break;
                    end
                end
                if found_predecessor
                    break;
                end
            end
        end
        return perm(sigma)
    end

    deck::Vector{PermGroupElem} = []
    degree_list = collect(1:degree(o))
    for i in degree_list
        candidate = candidate_for_deck(o, i)
        if (GapObj(candidate) != GAP.Globals.fail) && in_deck_group(o, candidate)
            push!(deck, candidate)
        end
    end
    return permutation_group(degree(o), deck)
end

function is_normal(o::Origami)
    return !(order(deck_group(o)) < degree(o)) 
end

# TODO implement AsNormalStoredOrigami