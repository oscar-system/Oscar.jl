
function orbit_weylgroup(L::LieAlgebraStructure, weight_vector_w::Vector{ZZRingElem})
  """
  operates weyl-group of type type and rank rank on vector weight_vector and returns list of vectors in orbit
  input and output weights in terms of w_i
  """
  # initialization
  weyl_group = GAP.Globals.WeylGroup(GAP.Globals.RootSystem(L.lie_algebra_gap))
  orbit_iterator = GAP.Globals.WeylOrbitIterator(weyl_group, GAP.Obj(Int.(weight_vector_w)))
  vertices = Vector{Int}[]

  # operate with the weylgroup on weight_vector
  GAPWrap.IsDoneIterator(orbit_iterator)
  while !(GAPWrap.IsDoneIterator(orbit_iterator))
    w = GAPWrap.NextIterator(orbit_iterator)
    push!(vertices, Vector{Int}(w))
  end

  return vertices
end

function get_dim_weightspace(
  L::LieAlgebraStructure, highest_weight::Vector{ZZRingElem}
)::Dict{Vector{ZZRingElem},Int}
  """
  Calculates dictionary with weights as keys and dimension of corresponding weightspace as value. GAP computes the 
  dimension for all positive weights. The dimension is constant on orbits of the weylgroup, and we can therefore 
  calculate the dimension of each weightspace. Returns weights in w_i
  """
  # calculate dimension for dominant weights with GAP
  root_system = GAP.Globals.RootSystem(L.lie_algebra_gap)
  dominant_char = GAP.Globals.DominantCharacter(root_system, GAP.Obj(Int.(highest_weight)))
  dominant_weights_w = map(weight -> ZZ.(weight), dominant_char[1])
  dominant_weights_dim = Int.(dominant_char[2])
  weightspaces = Dict{Vector{ZZRingElem},Int}()

  # calculate dimension for the rest by checking which positive weights lies in the orbit.
  for i in 1:length(dominant_weights_w)
    orbit_weights = orbit_weylgroup(L, dominant_weights_w[i])
    dim_weightspace = dominant_weights_dim[i]
    for weight in orbit_weights
      weightspaces[highest_weight - weight] = dim_weightspace
    end
  end
  return weightspaces
end

function convert_lattice_points_to_monomials(ZZx, lattice_points_weightspace)
  monomials = [
    finish(push_term!(MPolyBuildCtx(ZZx), ZZ(1), Int.(lattice_point))) for
    lattice_point in lattice_points_weightspace
  ]
  return monomials
end

function get_lattice_points_of_weightspace(
  weights_alpha::Vector{Vector{QQFieldElem}}, weight_alpha::Vector{QQFieldElem}
)
  """
  calculates all lattice points in a given weightspace for lie algebras
  input:
  weights: the operator weights in alpha_i
  weight: lambda - mu in alpha_i

  output: all lattice points with weight weight

  works by calculating all integer solutions to the following linear program:
  [   |              |    ]       [   x   ]      
  [weights[1]...weights[k]]   *   [   |   ]   =   weight 
  [   |              |    ]       [  res  ] 
  [   |              |    ]       [   |   ]
  where res[i] >= 0 for all i
  """
  n = length(weights_alpha)
  m = length(weight_alpha)
  A = zero_matrix(QQ, 2m + n, n)
  b = [zero(QQ) for _ in 1:(2m + n)]
  # equalities
  for i in 1:n
    w = matrix(QQ, m, 1, weights_alpha[i])
    A[1:m, i] = w
    A[(m + 1):(2m), i] = -w
  end
  b[1:m] = weight_alpha
  b[(m + 1):(2m)] = -weight_alpha
  # non-negativity
  for i in 1:n
    A[2m + i, i] = -1
  end

  return Vector{ZZRingElem}.(lattice_points(polyhedron(A, b)))
end
