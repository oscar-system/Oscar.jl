# Try to 'update' the base_ring of G
function map_entries(K::Field, G::MatrixGroup)
  g = dense_matrix_type(K)[]
  for h in gens(G)
    push!(g, map_entries(K, h.elm))
  end
  return matrix_group(g)
end

function is_reflection(g::MatrixGroupElem)
  return rank(g.elm - one(parent(g)).elm) == 1
end

function subgroup_of_reflections(G::MatrixGroup)
  g = elem_type(G)[]
  for c in conjugacy_classes(G)
    if is_reflection(representative(c))
      append!(g, collect(c))
    end
  end
  return matrix_group(degree(G), base_ring(G), g)
end

# Check if G contains reflections
function is_small_group(G::MatrixGroup)
  return !any(c -> is_reflection(representative(c)), conjugacy_classes(G))
end
