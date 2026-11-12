


@doc raw"""
    in_deck_group(o::Origami, sigma::PermGroupElem)

Checks whether or not the permutation sigma defines a deck transformation of the origami o.

# Examples
```jldoctest
julia> o = origami(cperm([1, 2, 3, 4, 5]), cperm())
Origami ((1,2,3,4,5),(), 5)

julia> in_deck_group(o, cperm(perm_group(o), [1, 2, 3, 4, 5]))
true

julia> in_deck_group(o, cperm(perm_group(o), [1, 2]))
false
```
"""
function in_deck_group(o::Origami, sigma::PermGroupElem)
  h = horizontal_perm(o)
  v = vertical_perm(o)
  return h == h^sigma && v == v^sigma
end

@doc raw"""
    deck_group(o::Origami)

Computes the group of deck transformations of the origami o as a covering of the once punctured torus. Note that the deck transformations can be seen as permutations of the squares of the origami.

# Examples
```jldoctest
julia> o = origami(cperm([1, 2, 3, 4, 5]), cperm())
Origami ((1,2,3,4,5),(), 5)

julia> deck_group(o)
Permutation group of degree 5 and order 5
```
"""
function deck_group(o::Origami)
  M = sub(perm_group(o), [horizontal_perm(o), vertical_perm(o)])[1]
  return centralizer(perm_group(o), M)[1]
end

function is_normal(o::Origami)
  # the deck group acts freely on each fibre, so its order is at most the degree
  return order(deck_group(o)) == degree(o)
end

# TODO implement AsNormalStoredOrigami
