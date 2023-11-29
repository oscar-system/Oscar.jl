
@doc raw"""
    orbit_weylgroup(L::LieAlgebraStructure, weight_w::Vector{ZZRingElem}) -> Vector{Vector{ZZRingElem}}

Computes the orbit of the weight `weight_w` given as coefficients to the fundamental weights
$\omega_i$ under the action of the Weyl group of the Lie algebra `L`.
"""
function orbit_weylgroup(L::LieAlgebraStructure, weight_w::Vector{ZZRingElem})
  # initialization
  weyl_group = GAP.Globals.WeylGroup(GAP.Globals.RootSystem(L.lie_algebra_gap))
  orbit_iterator = GAP.Globals.WeylOrbitIterator(weyl_group, GAP.Obj(Int.(weight_w)))
  vertices = Vector{ZZRingElem}[]

  # operate with the weylgroup on weight_vector
  while !(GAPWrap.IsDoneIterator(orbit_iterator))
    w = GAPWrap.NextIterator(orbit_iterator)
    push!(vertices, Vector{ZZRingElem}(w))
  end

  return vertices
end

@doc raw"""
    get_dim_weightspace(L::LieAlgebraStructure, highest_weight::Vector{ZZRingElem}) -> Dict{Vector{ZZRingElem},Int}

Computes the dimension of the weight spaces of the Lie algebra `L` module with highest weight `highest_weight`.
For all dominant weights, the dimension is computed with GAP. For the remaining weights, the dimension is
calculated by checking which dominant weight lies in the orbit of the weight under the action of the Weyl group.

The weights are given as coefficients to the fundamental weights $\omega_i$.
"""
function get_dim_weightspace(L::LieAlgebraStructure, highest_weight::Vector{ZZRingElem})
  # calculate dimension for dominant weights with GAP
  root_system = root_system_gap(L)
  dominant_char = GAP.Globals.DominantCharacter(root_system, GAP.Obj(Int.(highest_weight)))
  dominant_weights = map(weight -> ZZ.(weight), dominant_char[1])
  dominant_weights_dim = Int.(dominant_char[2])
  weightspaces = Dict{Vector{ZZRingElem},Int}()

  # calculate dimension for the rest by checking which dominant weight lies in the orbit
  for (dominant_weight, dim) in zip(dominant_weights, dominant_weights_dim)
    for weight in orbit_weylgroup(L, dominant_weight)
      weightspaces[highest_weight - weight] = dim
    end
  end
  return weightspaces
end

function convert_lattice_points_to_monomials(
  ZZx::ZZMPolyRing, lattice_points_weightspace::Vector{Vector{ZZRingElem}}
)
  return [ZZx([ZZ(1)], [lattice_point]) for lattice_point in lattice_points_weightspace]
end

@doc raw"""
    get_lattice_points_of_weightspace(weights::Vector{Vector{QQFieldElem}}, wght::Vector{QQFieldElem}) -> Vector{Vector{ZZRingElem}}

Calculates all lattice points in a given weightspace for a Lie algebra highest weight module.
This is equivalent to finding $\mathbb{Z}$-linear combinations of `weights` that equal `wght`.
All weights are given as coefficients to the simple roots $\alpha_i$.
"""
function get_lattice_points_of_weightspace(
  weights::Vector{Vector{QQFieldElem}}, wght::Vector{QQFieldElem}
)
  # calculate all integer solutions to the following linear program:
  # [   |              |    ]       [   x   ]      
  # [weights[1]...weights[k]]   *   [   |   ]   =   weight 
  # [   |              |    ]       [  res  ] 
  # [   |              |    ]       [   |   ]
  # where res[i] >= 0 for all i

  n = length(weights)
  m = length(wght)
  A = zero_matrix(QQ, 2m + n, n)
  b = [zero(QQ) for _ in 1:(2m + n)]
  # equalities
  for i in 1:n
    w = matrix(QQ, m, 1, weights[i])
    A[1:m, i] = w
    A[(m + 1):(2m), i] = -w
  end
  b[1:m] = wght
  b[(m + 1):(2m)] = -wght
  # non-negativity
  for i in 1:n
    A[2m + i, i] = -1
  end

  return Vector{ZZRingElem}.(lattice_points(polyhedron(A, b)))
end
