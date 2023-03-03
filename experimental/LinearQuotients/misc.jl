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

# Somehow one needs to look up powers of a given root of unity quite often...
# Assumes that zeta is a primitive l-th root of unity
function _powers_of_root_of_unity(zeta::FieldElem, l::Int)
  K = parent(zeta)
  powers_of_zeta = Dict{elem_type(K), Int}()
  t = one(K)
  for i = 0:l - 1
    powers_of_zeta[t] = i
    t *= zeta
  end
  @assert is_one(t) "Given element is not a $l-th root of unity"

  return powers_of_zeta
end

function is_subgroup_of_sl(G::MatrixGroup)
  for c in conjugacy_classes(G)
    if !is_one(det(representative(c).elm))
      return false
    end
  end
  return true
end
