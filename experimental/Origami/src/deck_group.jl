function in_deck_group(o::Origami, sigma::PermGroupElem)
  h = horizontal_perm(o)
  v = vertical_perm(o)
  return h == h^sigma && v == v^sigma
end

function deck_group(o::Origami)
  function candidate_for_deck(origami, j)
    h = horizontal_perm(origami)
    v = vertical_perm(origami)
    sym = perm_group(origami)
    sheets_to_visit = collect(2:degree(origami))
    sigma = fill(-1, degree(origami))
    sigma[1] = j
    while length(sheets_to_visit) > 0
      found_predecessor = false
      for i in sheets_to_visit
        for tau in (h, v)
          pos_iterator = searchsorted(sheets_to_visit, i^(tau^-1))
          if first(pos_iterator) > last(pos_iterator)
            found_predecessor = true
            sigma[i] = sigma[i^(tau^-1)]^tau
            deleteat!(sheets_to_visit, findfirst(item -> item == i,
              sheets_to_visit))
            break
          end
        end
        found_predecessor && break
      end
    end
    return perm(sym, sigma)
  end

  deck = PermGroupElem[]
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
