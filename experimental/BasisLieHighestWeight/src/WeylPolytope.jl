
function orbit_weylgroup(
  lie_algebra::LieAlgebraStructure, weight_vector_w::Vector{ZZRingElem}
)
  """
  operates weyl-group of type type and rank rank on vector weight_vector and returns list of vectors in orbit
  input and output weights in terms of w_i
  """
  # initialization
  weyl_group = GAP.Globals.WeylGroup(GAP.Globals.RootSystem(lie_algebra.lie_algebra_gap))
  weight_vector_w_int = convert(Vector{Int}, weight_vector_w)
  orbit_iterator = GAP.Globals.WeylOrbitIterator(weyl_group, GAP.Obj(weight_vector_w_int))
  vertices = []

  # operate with the weylgroup on weight_vector
  GAPWrap.IsDoneIterator(orbit_iterator)
  while !(GAPWrap.IsDoneIterator(orbit_iterator))
    w = GAPWrap.NextIterator(orbit_iterator)
    push!(vertices, Vector{Int}(w))
  end

  # return
  vertices = convert(Vector{Vector{Int}}, vertices)
  return vertices
end

function get_dim_weightspace(
  lie_algebra::LieAlgebraStructure, highest_weight::Vector{ZZRingElem}
)::Dict{Vector{ZZRingElem},Int}
  """
  Calculates dictionary with weights as keys and dimension of corresponding weightspace as value. GAP computes the 
  dimension for all positive weights. The dimension is constant on orbits of the weylgroup, and we can therefore 
  calculate the dimension of each weightspace. Returns weights in w_i
  """
  # calculate dimension for dominant weights with GAP
  root_system = GAP.Globals.RootSystem(lie_algebra.lie_algebra_gap)
  highest_weight_int = convert(Vector{Int}, highest_weight)
  result = GAP.Globals.DominantCharacter(root_system, GAP.Obj(highest_weight_int))
  dominant_weights_w = [map(Int, item) for item in result[1]]
  dominant_weights_dim = map(Int, result[2])
  dominant_weights_w = convert(Vector{Vector{ZZRingElem}}, dominant_weights_w)
  weightspaces = Dict{Vector{ZZRingElem},Int}()

  # calculate dimension for the rest by checking which positive weights lies in the orbit.
  for i in 1:length(dominant_weights_w)
    orbit_weights = orbit_weylgroup(lie_algebra, dominant_weights_w[i])
    dim_weightspace = dominant_weights_dim[i]
    for weight in orbit_weights
      weightspaces[highest_weight - weight] = dim_weightspace
    end
  end
  return weightspaces
end

function convert_lattice_points_to_monomials(ZZx, lattice_points_weightspace)
  #println("convert_lattice_points_to_monomials")
  #println("lattice_points_weightspace type: ", typeof(lattice_points_weightspace))
  #lat = [convert(Vector{Int}, convert(Vector{Int64}, lattice_point)) for lattice_point in lattice_points_weightspace]
  #println("before build")
  #monomials = [finish(push_term!(MPolyBuildCtx(ZZx), ZZ(1), p)) for p in lat]
  #          for lattice_point in lattice_points_weightspace]
  monomials = [
    finish(
      push_term!(
        MPolyBuildCtx(ZZx),
        ZZ(1),
        convert(Vector{Int}, convert(Vector{Int64}, lattice_point)),
      ),
    ) for lattice_point in lattice_points_weightspace
  ]
  #println("end convert_lattice_points_to_monomials")
  return monomials
end

function scale_weights_to_integers(
  weights_eps::Vector{Vector{QQFieldElem}}, weight_eps::Vector{QQFieldElem}
)
  # Extract all denominators from both structures
  denominators = [denominator(r) for w in weights_eps for r in w]
  append!(denominators, [denominator(r) for r in weight_eps])

  # Compute the LCM of all the denominators
  lcm_denominator = lcm(denominators...)

  # Scale the elements of weights_eps and weight_eps
  scaled_weights_eps = [[Int(lcm_denominator * r) for r in w] for w in weights_eps]
  scaled_weight_eps = [Int(lcm_denominator * r) for r in weight_eps]

  return scaled_weights_eps, scaled_weight_eps
end

function get_lattice_points_of_weightspace(
  weights_eps::Vector{Vector{QQFieldElem}},
  weight_eps::Vector{QQFieldElem},
  lie_type::Symbol,
)
  """
  calculates all lattice points in a given weightspace for a lie algebra of type type
  input:
  weights: the operator weights in eps_i
  weight: lambda - mu

  output: all lattice points with weight weight
  """
  # Rescale to integers
  scaled_weights_eps, scaled_weight_eps = scale_weights_to_integers(weights_eps, weight_eps)
  if lie_type in [:A, :G]
    return get_lattice_points_of_weightspace_A_G_n(scaled_weights_eps, scaled_weight_eps)
  else
    return get_lattice_points_of_weightspace_Xn(scaled_weights_eps, scaled_weight_eps)
  end
end

function get_lattice_points_of_weightspace_A_G_n(weights_eps, weight_eps)
  """
  calculates all monomials in a given weightspace for lie algebras that have type A or G
  input:
  weights: the operator weights in eps_i
  weight: lambda - mu

  output: all monomials with weight weight

  works by calculating all integer solutions to the following linear program:
  [    |              |     1 ]       [   |   ]      
  [weights[1]... weights[k] 1 ]   *   [  res  ]   =   weight 
  [    |              |     1 ]       [   |   ] 
  [    |              |     1 ]       [   x   ]
  where res[i] >= 0 for all i

  """
  n = length(weights_eps)
  m = length(weight_eps)
  A = zero_matrix(QQ, 2m + n, n + 1)
  b = [zero(QQ) for _ in 1:(2m + n)]
  # equalities
  for i in 1:n
    w = matrix(QQ, m, 1, weights_eps[i])
    A[1:m, i] = w
    A[(m + 1):(2m), i] = -w
  end
  A[1:m, n + 1] = [one(QQ) for _ in 1:m]
  A[(m + 1):(2m), n + 1] = [-one(QQ) for _ in 1:m]
  b[1:m] = weight_eps
  b[(m + 1):(2m)] = -weight_eps
  # non-negativity
  for i in 1:n
    A[2m + i, i] = -1
  end

  return [point[1:n] for point in lattice_points(polyhedron(A, b))]
end

function get_lattice_points_of_weightspace_Xn(weights_eps, weight_eps)
  """
  calculates all lattice points in a given weightspace for lie algebras that don't have type A or G
  input:
  weights: the operator weights in eps_i
  weight: lambda - mu

  output: all lattice points with weight weight

  works by calculating all integer solutions to the following linear program:
  [   |              |    ]       [   x   ]      
  [weights[1]...weights[k]]   *   [   |   ]   =   weight 
  [   |              |    ]       [  res  ] 
  [   |              |    ]       [   |   ]
  where res[i] >= 0 for all i
  """
  n = length(weights_eps)
  m = length(weight_eps)
  A = zero_matrix(QQ, 2m + n, n)
  b = [zero(QQ) for _ in 1:(2m + n)]
  # equalities
  for i in 1:n
    w = matrix(QQ, m, 1, weights_eps[i])
    A[1:m, i] = w
    A[(m + 1):(2m), i] = -w
  end
  b[1:m] = weight_eps
  b[(m + 1):(2m)] = -weight_eps
  # non-negativity
  for i in 1:n
    A[2m + i, i] = -1
  end

  return lattice_points(polyhedron(A, b))
end
