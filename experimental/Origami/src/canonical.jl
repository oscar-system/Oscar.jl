"""
Return a canonical, immutable key for a connected pair of permutations.

The permutations `h` and `v` use Julia labels 1:n. Two connected pairs have
the same key exactly when they are simultaneously conjugate.

The result contains first the canonical horizontal permutation and then the
canonical vertical permutation.
"""
function canonical_origami_key(o::Origami)
  n = degree(o)
  h = Vector(horizontal_perm(o))
  v = Vector(vertical_perm(o))

  # Reused workspace.
  old_to_new = zeros(Int, n)
  queue = Vector{Int}(undef, n)
  candidate_h = Vector{Int}(undef, n)
  candidate_v = Vector{Int}(undef, n)

  # Compute the canonical candidate associated with one root.
  # The 'root' is the arbitrary square (with the original labeling of o) which will
  # be 1 in the new labeling. Outgoing from this root, the rest of the labeling is
  # built in a deterministic and thus unique way.
  function candidate_for_root!(root::Int)
    fill!(old_to_new, 0)

    queue[1] = root
    old_to_new[root] = 1

    head = 1
    tail = 1

    while head <= tail
      x = queue[head]
      head += 1

      y = h[x] # move along horizontally
      if old_to_new[y] == 0 # if the square y (to the right of x) is not labeled yet
        tail += 1 # increase the counter
        queue[tail] = y # add y to our queue
        old_to_new[y] = tail # fix the new label for y
      end

      y = v[x] # move along vertically
      if old_to_new[y] == 0 # if the square y (above x) is not labeled yet
        tail += 1 # increase the counter
        queue[tail] = y # add y to our queue
        old_to_new[y] = tail # fix the new label for y
      end
    end

    # build the new candidate permutation
    for old_x in 1:n
      new_x = old_to_new[old_x]
      candidate_h[new_x] = old_to_new[h[old_x]]
      candidate_v[new_x] = old_to_new[v[old_x]]
    end
    return
  end

  # Use root 1 to initialize concrete best arrays.
  candidate_for_root!(1)

  best_h = copy(candidate_h)
  best_v = copy(candidate_v)

  for root in 2:n
    candidate_for_root!(root)

    # choose the lexicographically smallest tuple of permutations
    # (the lexicographic ordering applied to the list representations)
    if (candidate_h, candidate_v) < (best_h, best_v)
      copyto!(best_h, candidate_h)
      copyto!(best_v, candidate_v)
    end
  end

  return CanonicalOrigamiKey(best_h, best_v)
end

function normal_form(o::Origami)
  k = canonical_origami_key(o)
  h = perm(k.h)
  v = perm(k.v)
  return origami_disconnected(h, v, degree(o))
end
