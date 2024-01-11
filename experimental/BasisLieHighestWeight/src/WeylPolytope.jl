
@doc raw"""
    orbit_weylgroup(L::LieAlgebraStructure, weight_w::Vector{ZZRingElem}) -> Vector{Vector{ZZRingElem}}

Computes the orbit of the weight `weight_w` given as coefficients to the fundamental weights
$\omega_i$ under the action of the Weyl group of the Lie algebra `L`.
"""
function orbit_weylgroup(L::LieAlgebraStructure, weight_w::Vector{ZZRingElem})
  # initialization
  weyl_group = GAPWrap.WeylGroup(GAPWrap.RootSystem(L.lie_algebra_gap))
  orbit_iterator = GAPWrap.WeylOrbitIterator(weyl_group, GAP.Obj(Int.(weight_w)))
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
    get_lattice_points_of_weightspace(weight_roots::Vector{Vector{QQFieldElem}}, L::LieAlgebraStructure, weight::Vector{ZZRingElem}, highest_weight::Vector{ZZRingElem},


Calculates all lattice points in a given weightspace for a Lie algebra highest weight module.
This is equivalent to finding $\mathbb{Z}$-linear combinations of `weights` that equal `wght`.
All weights are given as coefficients to the simple roots $\alpha_i$.
"""
function get_lattice_points_of_weightspace(
  weight_roots::Vector{Vector{QQFieldElem}}, 
  L::LieAlgebraStructure, 
  weight::Vector{ZZRingElem}, 
  highest_weight::Vector{ZZRingElem}, 
  zero_coordinates::Set{Any}
)
  # calculate all integer solutions to the following linear program:
  # [   |              |    ]       [   x   ]      
  # [weight_roots[1]...weight_roots[k]]   *   [   |   ]   =   weight
  # [   |              |    ]       [  res  ] 
  # [   |              |    ]       [   |   ]
  # where res[i] >= 0 for all i
  wght = w_to_alpha(L, weight )
  n = length(weight_roots)
  m = length(wght)
 
  #zero_coordinates = 0
  #zero_coordinates = compute_zero_coordinates(weight_roots, L, weight, highest_weight) 
  A = zero_matrix(QQ, 2m + n + length(zero_coordinates), n)
  b = [zero(QQ) for _ in 1:(2m + n + length(zero_coordinates))]

  # equalities
  for i in 1:n
    w = matrix(QQ, m, 1, weight_roots[i])
    A[1:m, i] = w
    A[(m + 1):(2m), i] = -w
  end

  b[1:m] = wght
  b[(m + 1):(2m)] = -wght
  # non-negativity
  for i in 1:n
    A[2m + i, i] = -1
  end
  for (j, i) in enumerate(zero_coordinates) 
    A[2m + n + j, i] = 1
    b[2m + n + j] = 0
  end
  sol = Vector{ZZRingElem}.(lattice_points(polyhedron(A, b)))
  return sol
end

@doc raw"""
compute_zero_coordinates(
  weight_roots::Vector{Vector{QQFieldElem}}, L::LieAlgebraStructure, weight::Vector{ZZRingElem}, highest_weight::Vector{ZZRingElem},
)
Uses the action on the generator to determine coordinates which always act by 0 
on the generator, these are the root vectors which are supported on a complement of the support of the highest weight. 
This has to be read from right to left, each possible action with a root vector increases the support of the so far generated module.
 This might speed up the function get_lattice_points_of_weightspace if the roots to the right end of the birational sequence are small.
"""

function compute_zero_coordinates(
  weight_roots::Vector{Vector{QQFieldElem}}, L::LieAlgebraStructure, weight::Vector{ZZRingElem}, highest_weight::Vector{ZZRingElem},
)
  n = length(weight_roots)
  m = length(weight)
  a = cartan_matrix(L)
  non_zeros = Set()
  all = Set()
  for i in 1:m 
    push!(all, i)
  end
  non_zeros = Set(findall(!iszero , highest_weight))

  zero_coordinates = Set()
  c = n
  check = false
  collect =Set()
  while c > 0 && !issubset(all, non_zeros)
    for i in 1:m
      if (weight_roots[c][i] != 0) && ( i in non_zeros)
        check = true
        push!(collect, i)
      end
    end  
    for i in collect
      for j in 1:m
        if (a[i,j] != 0)
          push!(non_zeros, j)
        end
      end   
    end   
    if (!check)
      push!(zero_coordinates, c)
    end
    c = c-1
    check = false
  end 
  return zero_coordinates
end
