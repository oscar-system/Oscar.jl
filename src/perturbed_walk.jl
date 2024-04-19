###############################################################
# Perturbed-version of the Groebner Walk.
###############################################################

function perturbed_walk(
  G::Oscar.IdealGens,
  start::MonomialOrdering,
  target::MonomialOrdering,
  p::Int
)
  @vprintln :groebner_walk "perturbed_walk results"
  @vprintln :groebner_walk "Crossed Cones in: "

  R = base_ring(G)
  S = canonical_matrix(start)
  T = canonical_matrix(target)

  current_weight = perturbed_vector(G, S, p)

  while !same_cone(G, target)
    target_weight = perturbed_vector(G, T, p)
    next_target = matrix_ordering(R, add_weight_vector(target_weight, T))
    G = standard_walk(G, next_target, current_weight, target_weight)

    p = p - 1
    current_weight = target_weight
    S = next_target
  end

  return G
end

perturbed_walk(G::Oscar.IdealGens, S::Matrix{Int}, T::Matrix{Int}, p::Int) = perturbed_walk(G, matrix_ordering(base_ring(G), S), matrix_ordering(base_ring(G), T), p)

# computes a p-perturbed vector from the matrix M.
function perturbed_vector(G::Oscar.IdealGens, M::ZZMatrix, p::Int)
  n = size(M, 1)
  rows = [M[i, :] for i in 1:p]

  m = maximum.(Ref(abs), rows)
  m_sum = sum(m)
  max_deg = maximum(total_degree.(G)) # TODO: I think this is total degree

  e = max_deg * m_sum + 1
  w = M[1, :] * e^(p - 1)
  for i in 2:p
    w += e^(p - i) * M[i, :]
  end

  w = sum(rows .* (p .- Vector(1:p)))

  return convert_bounding_vector(w)
end