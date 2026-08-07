function normal_form(o::Origami)
  n = degree(o)

  sym = perm_group(o)
  x = sym(horizontal_perm(o))
  y = sym(vertical_perm(o))

  # Find points which minimize the lengths of the cycles in which they occur.
  # This can greatly reduce the number of breadths-first searches below.
  minimal_cycle_lengths = (n, n)
  minimize_cycle_lengths = Int[]
  for i in 1:n
    cycle_lengths = (cycle_length(x, i), cycle_length(y, i))
    if cycle_lengths == minimal_cycle_lengths
      push!(minimize_cycle_lengths, i)
    elseif cycle_lengths < minimal_cycle_lengths
      empty!(minimize_cycle_lengths)
      push!(minimize_cycle_lengths, i)
      minimal_cycle_lengths = cycle_lengths
    end
  end

  G = PermGroupElem[]

  # TODO this can be completely refactored. this is the exact same as
  # normalform_conjugators apart from looping over different lists
  L = fill(0, n)
  Q = Int[]
  for i in minimize_cycle_lengths
    fill!(L, 0)
    empty!(Q)
    push!(Q, i)
    numSeen = 1
    L[i] = 1
    while numSeen < n
      v = popfirst!(Q)
      wx = v^x
      wy = v^y
      if L[wx] == 0
        push!(Q, wx)
        numSeen += 1
        L[wx] = numSeen
      end
      if L[wy] == 0
        push!(Q, wy)
        numSeen += 1
        L[wy] = numSeen
      end
    end
    push!(G, perm(sym, L))
  end

  G2 = [(x^i, y^i) for i in G]
  # no need to check if surface connected
  min_entry = minimum(G2)
  return origami_disconnected(min_entry[1], min_entry[2], n)
end

"""
Return whether `(h1, v1)` is lexicographically smaller than`(h2, v2)`.
"""
function permutation_pair_isless(h1::Vector{Int}, v1::Vector{Int}, h2::Vector{Int},v2::Vector{Int})::Bool

  for i in eachindex(h1, h2)
    h1[i] < h2[i] && return true
    h1[i] > h2[i] && return false
  end

  for i in eachindex(v1, v2)
    v1[i] < v2[i] && return true
    v1[i] > v2[i] && return false
  end

  return false
end

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
  function candidate_for_root!(root::Int)
    fill!(old_to_new, 0)

    queue[1] = root
    old_to_new[root] = 1

    head = 1
    tail = 1

    while head <= tail
      x = queue[head]
      head += 1

      y = h[x]
      if old_to_new[y] == 0
        tail += 1
        queue[tail] = y
        old_to_new[y] = tail
      end

      y = v[x]
      if old_to_new[y] == 0
        tail += 1
        queue[tail] = y
        old_to_new[y] = tail
      end
    end

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

    if permutation_pair_isless(candidate_h, candidate_v, best_h, best_v)
      copyto!(best_h, candidate_h)
      copyto!(best_v, candidate_v)
    end
  end

  return CanonicalOrigamiKey(best_h, best_v)
end
