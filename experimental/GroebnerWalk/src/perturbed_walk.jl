@doc raw"""
    perturbed_walk(G::Oscar.IdealGens, start::MonomialOrdering, target::MonomialOrdering, p::Int)

Compute a reduced Groebner basis w.r.t. to a monomial order by converting it using 
the Groebner Walk with the algorithm proposed by [Tra00](@cite).

# Arguments
- `G::Oscar.IdealGens`: Groebner basis of an ideal with respect to a starting monomial order.
- `target::MonomialOrdering`: monomial order one wants to compute a Groebner basis for.
- `start::MonomialOrdering`: monomial order to begin the conversion.
- `p::Int`: degree of perturbation
"""
function perturbed_walk(
  G::Oscar.IdealGens,
  start::MonomialOrdering,
  target::MonomialOrdering,
  p::Int
)
  @vprintln :groebner_walk "Results for perturbed_walk"
  @vprintln :groebner_walk "Crossed Cones in: "

  R = base_ring(G)
  S = canonical_matrix(start)
  T = canonical_matrix(target)

  current_weight = perturbed_vector(G, S, p)

  while !same_cone(G, target)
    # @v_do :groebner_walk steps += 1
    @vprintln :groebner_walk current_weight
    @vprintln :groebner_walk 2 "Next matrix for monomial order:"
    @vprintln :groebner_walk 2 S

    target_weight = perturbed_vector(G, T, p)
    next_target = weight_ordering(target_weight, target)
    G = standard_walk(Oscar.IdealGens, G, next_target, current_weight, target_weight)

    @vprintln :groebner_walk 2 G
    @vprintln :groebner_walk "======="

    p = p - 1
    current_weight = target_weight
  end

  return gens(G)
end

perturbed_walk(
  G::Oscar.IdealGens,
  S::ZZMatrix,
  T::ZZMatrix,
  p::Int
 ) = perturbed_walk(G, matrix_ordering(base_ring(G), S), matrix_ordering(base_ring(G), T), p)

@doc raw"""
    perturbed_vector(G::Oscar.IdealGens, M::ZZMatrix, p::Int)

Compute a perturbed vector using a matrix `M` representing some monomial order
for one iteration of the Groebner walk according to [Tra00](@cite), Thm. 3.1.

# Arguments
- `G::Oscar.IdealGens`: Groebner basis of an ideal with respect to a starting monomial order.
- `M::ZZMatrix`: matrix representing a monomial order
- `p::Int`: number of rows of `M` to use for the perturbation
"""
function perturbed_vector(G::Oscar.IdealGens, M::ZZMatrix, p::Int)
  n = size(M, 1)
  rows = [M[i, :] for i in 1:p]

  # Calculate upper degree bound for the target Groebner basis
  m = maximum.(Ref(abs), rows) |> maximum
  max_deg = maximum(total_degree.(G)) 
  e = max_deg * m + 1
  @vprint :groebner_walk 5 "Upper degree bound: "
  @vprintln :groebner_walk 5 e

  #w = sum(row * e^(p-i) for (i,row) in enumerate(rows))
  w = zeros(ZZRingElem, n)
  d = 1
  for i in 1:p
    w += row[end+1-i] * d
    d *= e
  end

  @vprint :groebner_walk 3 "Perturbed vector: "
  @vprintln :groebner_walk 3 w

  return convert_bounding_vector(w)
end

