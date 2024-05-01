#### PARAMETRISATION IN PROBABLITY COORDINATES #### 

function monomial_parametrization(pm::PhylogeneticModel, states::Dict{Int, Int})
  gr = graph(pm)
  tr_mat = transition_matrices(pm)
  root_dist = root_distribution(pm)

  r = root(gr)
  monomial = root_dist[states[r]]
  for edge in edges(gr)
    stateParent = states[src(edge)]
    stateChild = states[dst(edge)]
    monomial = monomial * tr_mat[edge][stateParent, stateChild]
  end

  return monomial
end
function monomial_parametrization(pm::GroupBasedPhylogeneticModel, states::Dict{Int, Int})
  monomial_parametrization(phylogenetic_model(pm), states)
end

function probability_parametrization(pm::PhylogeneticModel, leaves_states::Vector{Int})
  gr = graph(pm)
  int_nodes = interior_nodes(gr)
  lvs_nodes = leaves(gr)
  n_states = number_states(pm)

  interior_indices = collect.(Iterators.product([collect(1:n_states) for _ in int_nodes]...))  
  nodes_states = Dict(lvs_nodes[i] => leaves_states[i] for i in 1:length(lvs_nodes))

  poly = 0
  # Might be useful in the future to use a polynomial ring context
  for labels in interior_indices
    for (int_node, label) in zip(int_nodes, labels)
      nodes_states[int_node] = label
    end
    poly = poly + monomial_parametrization(pm, nodes_states)
  end 
  return poly
end 
function probability_parametrization(pm::GroupBasedPhylogeneticModel, leaves_states::Vector{Int})
  probability_parametrization(phylogenetic_model(pm), leaves_states)
end
 
@doc raw"""
    probability_map(pm::PhylogeneticModel)    

Create a parametrization for a `PhylogeneticModel` of type `Dictionary`.
Iterates through all possible states of the leaf random variables and calculates their corresponding probabilities using the root distribution and laws of conditional independence. Returns a dictionary of polynomials indexed by the states. Uses auxiliary function `monomial_parametrization(pm::PhylogeneticModel, states::Dict{Int, Int})` and `probability_parametrization(pm::PhylogeneticModel, leaves_states::Vector{Int})`. 

# Examples
```jldoctest
julia>
```
"""
function probability_map(pm::PhylogeneticModel)
  lvs_nodes = leaves(graph(pm))
  n_states = number_states(pm)

  leaves_indices = collect.(Iterators.product([collect(1:n_states) for _ in lvs_nodes]...))
  probability_coordinates = Dict(Tuple(leaves_states) => probability_parametrization(pm, leaves_states) for leaves_states in leaves_indices)
  return probability_coordinates
end
function probability_map(pm::GroupBasedPhylogeneticModel)
  probability_map(phylogenetic_model(pm))
end


#### FOURIER PARAMETRISATION #### 
  
function monomial_fourier(pm::GroupBasedPhylogeneticModel, leaves_states::Vector{Int})
  gr = graph(pm)
  param = fourier_parameters(pm)
  monomial = 1
  for edge in edges(gr)
    dsc = vertex_descendants(dst(edge), gr, [])
    elem = group_sum(pm, leaves_states[dsc])
    monomial = monomial * param[edge][which_group_element(pm, elem)]
  end
  return monomial
end
  
function fourier_parametrization(pm::GroupBasedPhylogeneticModel, leaves_states::Vector{Int})
  S = fourier_ring(pm)
  if is_zero_group_sum(pm, leaves_states) 
    poly = monomial_fourier(pm, leaves_states)
  else 
    poly = S(0)
  end

  return poly
end 


@doc raw"""
    fourier(pm::GroupBasedPhylogeneticModel)    

Create a parametrization for a `GroupBasedPhylogeneticModel` of type `Dictionary`.
Iterates through all possible states of the leaf random variables and calculates their corresponding probabilities using group actions and laws of conditional independence. Returns a dictionary of polynomials indexed by the states. Uses auxiliary function `monomial_fourier(pm::GroupBasedPhylogeneticModel, leaves_states::Vector{Int})` and `fourier_parametrization(pm::GroupBasedPhylogeneticModel, leaves_states::Vector{Int})`. 

# Examples
```jldoctest
julia>
```
"""
function fourier_map(pm::GroupBasedPhylogeneticModel)
  lvs_nodes = leaves(graph(pm))
  n_states = number_states(pm)

  leaves_indices = collect.(Iterators.product([collect(1:n_states) for _ in lvs_nodes]...))
  fourier_coordinates = Dict(Tuple(leaves_states) => fourier_parametrization(pm, leaves_states) for leaves_states in leaves_indices)
  return fourier_coordinates
end


#### SPECIALIZED FOURIER TRANSFORM MATRIX ####


@doc raw"""
    specialized_fourier_transform(pm::GroupBasedPhylogeneticModel, p_equivclasses::Dict{Vector{Vector{Int64}}, QQMPolyRingElem}, f_equivclasses::Dict{Vector{Vector{Int64}}, QQMPolyRingElem})    

Reparametrize between a model specification in terms of probability and Fourier cooordinates.

# Examples
```jldoctest
julia>
```
"""
function specialized_fourier_transform(pm::GroupBasedPhylogeneticModel, p_equivclasses::Dict{Vector{Vector{Int64}}, QQMPolyRingElem}, f_equivclasses::Dict{Vector{Vector{Int64}}, QQMPolyRingElem})
  R = probability_ring(pm)
  ns = number_states(pm)

  np = length(p_equivclasses)
  nq = length(f_equivclasses) - 1

  ## We need to sort the equivalence classes: both inside each class as well as the collection of classes. 
  p_equivclasses_sorted = collect(keys(p_equivclasses))
  [sort!(p_eqclass) for p_eqclass in p_equivclasses_sorted]
  sort!(p_equivclasses_sorted)

  f_equivclasses_sorted = collect(keys(filter(x -> !is_zero(x.second), f_equivclasses)))
  [sort!(f_eqclass) for f_eqclass in f_equivclasses_sorted]
  sort!(f_equivclasses_sorted)

  H = R.(hadamard(matrix_space(ZZ, ns, ns)))

  specialized_ft_matrix = R.(Int.(zeros(nq, np)))
  for i in 1:nq
    current_fourier_classes = f_equivclasses_sorted[i]
    for j in 1:np
      current_prob_classes = p_equivclasses_sorted[j]
      current_entriesin_M = [prod([H[y,x] for (x,y) in zip(p,q)]) for p in current_prob_classes, q in current_fourier_classes]
      specialized_ft_matrix[i,j] = R.(1//(length(current_prob_classes)*length(current_fourier_classes))*sum(current_entriesin_M))
    end
  end
  return specialized_ft_matrix
end

@doc raw"""
    inverse_specialized_fourier_transform(pm::GroupBasedPhylogeneticModel, p_equivclasses::Dict{Vector{Vector{Int64}}, QQMPolyRingElem},f_equivclasses::Dict{Vector{Vector{Int64}}, QQMPolyRingElem})    

Reparametrize between a model specification in terms of Fourier and probability cooordinates.

# Examples
```jldoctest
julia>
```
"""
function inverse_specialized_fourier_transform(pm::GroupBasedPhylogeneticModel, p_equivclasses::Dict{Vector{Vector{Int64}}, QQMPolyRingElem},f_equivclasses::Dict{Vector{Vector{Int64}}, QQMPolyRingElem})
  R = probability_ring(pm)
  ns = number_states(pm)

  np = length(p_equivclasses)
  nq = length(f_equivclasses) - 1

  ## We need to sort the equivalence classes: both inside each class as well as the collection of classes. 
  p_equivclasses_sorted = collect(keys(p_equivclasses))
  [sort!(p_eqclass) for p_eqclass in p_equivclasses_sorted]
  sort!(p_equivclasses_sorted)

  f_equivclasses_sorted = collect(keys(filter(x -> !is_zero(x.second), f_equivclasses)))
  [sort!(f_eqclass) for f_eqclass in f_equivclasses_sorted]
  sort!(f_equivclasses_sorted)

  H = R.(hadamard(matrix_space(ZZ, ns, ns)))
  Hinv = 1//ns * H 

  inverse_spec_ft_matrix = R.(Int.(zeros(np, nq)))
  for i in 1:np
    current_prob_class = p_equivclasses_sorted[i]
    for j in 1:nq
      current_fourier_class = f_equivclasses_sorted[j]
      current_entriesin_Minv = [prod([Hinv[x,y] for (x,y) in zip(p,q)]) for p in current_prob_class, q in current_fourier_class] 
      inverse_spec_ft_matrix[i,j] = R.(sum(current_entriesin_Minv))
    end
  end
  return inverse_spec_ft_matrix
end
