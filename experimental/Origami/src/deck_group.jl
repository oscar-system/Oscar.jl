function in_deck_group(o::Origami, sigma::PermGroupElem)
  h = horizontal_perm(o)
  v = vertical_perm(o)
  for i in 1:degree(o)
    if (i^sigma)^h != (i^h)^sigma || (i^sigma)^v != (i^v)^sigma
      return false
    end
  end
  return true
end

function deck_group(o::Origami)
  function candidate_for_deck(origami, j)
    sheets_to_visit = 2:degree(origami)
    sigma = [j]
    append!(sigma, fill(-1, degree(origami) - 1))
    while length(sheets_to_visit) > 0
      found_predecessor = false
      h = horizontal_perm(origami)
      v = vertical_perm(origami)
      for i in sheets_to_visit
        for tau in [h, v]
          pos_iterator = searchsorted(sheets_to_visit, i^(tau^-1))
          if first(pos_iterator) > last(pos_iterator)
            found_predecessor = true
            sigma[i] = sigma[i ^ (tau ^ -1)]^tau
            deleteat!(sheets_to_visit, findfirst(item -> item == i,
              sheets_to_visit))
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
  for i in 1:degree(o)
    candidate = candidate_for_deck(o, i)
    if GapObj(candidate) != GAP.Globals.fail && in_deck_group(o, candidate)
      push!(deck, candidate)
    end
  end
  return permutation_group(degree(o), deck)
end

function is_normal(o::Origami)
  return !(order(deck_group(o)) < degree(o))
end

# TODO implement AsNormalStoredOrigami
