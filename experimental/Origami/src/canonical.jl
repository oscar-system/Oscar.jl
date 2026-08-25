@doc raw"""
      normal_form(o::Origami)

Computes a normal form of a given origami. This normal form has the property that two origamis are equivalent if and only if their normal forms are equal. Note that this method does not copy any previously computed attributes from the old origami to the new one.

# Examples
```jldoctest
julia> o = origami(cperm([1,3,2,4]),cperm([2,3]))
Origami ((1,3,2,4),(2,3), 4)

julia> normal_form(o)
Origami ((1,2,3,4),(3,4), 4)
```
"""
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
