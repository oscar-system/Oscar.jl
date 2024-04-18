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
  
  perturbed_walk(G::Oscar.IdealGens, S::Matrix{Int}, T::Matrix{Int}, p::Int) = perturbed_walk(G, matrix_ordering(base_ring(G),S), matrix_ordering(base_ring(G),T), p)
  