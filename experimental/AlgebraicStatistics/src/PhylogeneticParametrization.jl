####################################################
#### PARAMETRIZATION IN PROBABILITY COORDINATES ####
####################################################

function monomial_parametrization(pm::PhylogeneticModel, states::Dict{Int, Int})
  gr = graph(pm)
  par_gens = parameter_ring_gens(pm)
  tr_mat = trans_matrix(pm)
  root_dist = root_distribution(pm)

  r = root(gr)
  monomial = root_dist[states[r]]
  
  for (i, edge) in enumerate(edges(gr))
    state_parent = states[src(edge)]
    state_child = states[dst(edge)]
    # get the symbolfrom the transition matrix signature
    sym = tr_mat[state_parent, state_child]
    # we can try and avoid this here if this becomes a bottle neck
    # it's related to the comment below about using a polynomial context
    # i.e., this is just the adding of exponents which are integer vectors
    # which would be much faster than polynomial multiplication
    monomial = monomial * par_gens[(sym, i)]
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

Iterate through all possible states of the leaf random variables and calculates their corresponding probabilities using the root distribution and laws of conditional independence. Return a dictionary of polynomials indexed by the states. Use auxiliary function `monomial_parametrization(pm::PhylogeneticModel, states::Dict{Int, Int})` and `probability_parametrization(pm::PhylogeneticModel, leaves_states::Vector{Int})`. 

# Examples
```jldoctest
julia> pm = jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]));

julia> p = probability_map(pm)
Dict{Tuple{Vararg{Int64}}, QQMPolyRingElem} with 64 entries:
  (1, 2, 1) => 1//4*a[1]*a[3]*b[2] + 1//4*a[2]*b[1]*b[3] + 1//2*b[1]*b[2]*b[3]
  (3, 1, 1) => 1//4*a[1]*b[2]*b[3] + 1//4*a[2]*a[3]*b[1] + 1//2*b[1]*b[2]*b[3]
  (4, 4, 2) => 1//4*a[1]*a[2]*b[3] + 1//4*a[3]*b[1]*b[2] + 1//2*b[1]*b[2]*b[3]
  (1, 2, 3) => 1//4*a[1]*b[2]*b[3] + 1//4*a[2]*b[1]*b[3] + 1//4*a[3]*b[1]*b[2] …
  (3, 1, 3) => 1//4*a[1]*a[3]*b[2] + 1//4*a[2]*b[1]*b[3] + 1//2*b[1]*b[2]*b[3]
  (3, 2, 4) => 1//4*a[1]*b[2]*b[3] + 1//4*a[2]*b[1]*b[3] + 1//4*a[3]*b[1]*b[2] …
  (3, 2, 1) => 1//4*a[1]*b[2]*b[3] + 1//4*a[2]*b[1]*b[3] + 1//4*a[3]*b[1]*b[2] …
  (2, 1, 4) => 1//4*a[1]*b[2]*b[3] + 1//4*a[2]*b[1]*b[3] + 1//4*a[3]*b[1]*b[2] …
  (3, 2, 3) => 1//4*a[1]*a[3]*b[2] + 1//4*a[2]*b[1]*b[3] + 1//2*b[1]*b[2]*b[3]
  (2, 1, 1) => 1//4*a[1]*b[2]*b[3] + 1//4*a[2]*a[3]*b[1] + 1//2*b[1]*b[2]*b[3]
  (1, 3, 2) => 1//4*a[1]*b[2]*b[3] + 1//4*a[2]*b[1]*b[3] + 1//4*a[3]*b[1]*b[2] …
  (1, 4, 2) => 1//4*a[1]*b[2]*b[3] + 1//4*a[2]*b[1]*b[3] + 1//4*a[3]*b[1]*b[2] …
  (2, 1, 3) => 1//4*a[1]*b[2]*b[3] + 1//4*a[2]*b[1]*b[3] + 1//4*a[3]*b[1]*b[2] …
  (2, 2, 4) => 1//4*a[1]*a[2]*b[3] + 1//4*a[3]*b[1]*b[2] + 1//2*b[1]*b[2]*b[3]
  (4, 3, 4) => 1//4*a[1]*a[3]*b[2] + 1//4*a[2]*b[1]*b[3] + 1//2*b[1]*b[2]*b[3]
  (2, 2, 1) => 1//4*a[1]*a[2]*b[3] + 1//4*a[3]*b[1]*b[2] + 1//2*b[1]*b[2]*b[3]
  (4, 4, 4) => 1//4*a[1]*a[2]*a[3] + 3//4*b[1]*b[2]*b[3]
  (4, 3, 1) => 1//4*a[1]*b[2]*b[3] + 1//4*a[2]*b[1]*b[3] + 1//4*a[3]*b[1]*b[2] …
  (3, 3, 2) => 1//4*a[1]*a[2]*b[3] + 1//4*a[3]*b[1]*b[2] + 1//2*b[1]*b[2]*b[3]
  ⋮         => ⋮
```
"""
function probability_map(pm::PhylogeneticModel)
  lvs_nodes = leaves(graph(pm))
  n_states = number_states(pm)

  leaves_indices = collect.(Iterators.product([collect(1:n_states) for _ in lvs_nodes]...))
  probability_coordinates = Dict{Tuple{Vararg{Int64}}, QQMPolyRingElem}(Tuple(leaves_states) => probability_parametrization(pm, leaves_states) for leaves_states in leaves_indices)
  return probability_coordinates
end

function probability_map(pm::GroupBasedPhylogeneticModel)
  probability_map(phylogenetic_model(pm))
end


################################################
#### PARAMETRIZATION IN FOURIER COORDINATES ####
################################################
  
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
    fourier_map(pm::GroupBasedPhylogeneticModel)    

Create a parametrization for a `GroupBasedPhylogeneticModel` of type `Dictionary`.

Iterate through all possible states of the leaf random variables and calculates their corresponding probabilities using group actions and laws of conditional independence. Return a dictionary of polynomials indexed by the states. Use auxiliary function `monomial_fourier(pm::GroupBasedPhylogeneticModel, leaves_states::Vector{Int})` and `fourier_parametrization(pm::GroupBasedPhylogeneticModel, leaves_states::Vector{Int})`. 

# Examples
```jldoctest
julia> pm = jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]));

julia> q = fourier_map(pm)
Dict{Tuple{Vararg{Int64}}, QQMPolyRingElem} with 64 entries:
  (1, 2, 1) => 0
  (3, 1, 1) => 0
  (4, 4, 2) => 0
  (1, 2, 3) => 0
  (3, 1, 3) => x[2, 1]*x[1, 2]*x[3, 2]
  (3, 2, 4) => x[1, 2]*x[2, 2]*x[3, 2]
  (3, 2, 1) => 0
  (2, 1, 4) => 0
  (3, 2, 3) => 0
  (2, 1, 1) => 0
  (1, 3, 2) => 0
  (1, 4, 2) => 0
  (2, 1, 3) => 0
  (2, 2, 4) => 0
  (4, 3, 4) => 0
  (2, 2, 1) => x[3, 1]*x[1, 2]*x[2, 2]
  (4, 4, 4) => 0
  (4, 3, 1) => 0
  (3, 3, 2) => 0
  ⋮         => ⋮
```
"""
function fourier_map(pm::GroupBasedPhylogeneticModel)
  lvs_nodes = leaves(graph(pm))
  n_states = number_states(pm)

  leaves_indices = collect.(Iterators.product([collect(1:n_states) for _ in lvs_nodes]...))
  fourier_coordinates = Dict{Tuple{Vararg{Int64}}, QQMPolyRingElem}(Tuple(leaves_states) => fourier_parametrization(pm, leaves_states) for leaves_states in leaves_indices)
  return fourier_coordinates
end


#####################################
#### COMPUTE EQUIVALENCE CLASSES ####
#####################################

@doc raw"""
    compute_equivalent_classes(parametrization::Dict{Tuple{Vararg{Int64}}, QQMPolyRingElem})

Given the parametrization of a `PhylogeneticModel`, cancel all duplicate entries and return equivalence classes of states which are attached the same probabilities.

# Examples
```jldoctest
julia> pm = jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]));

julia> p = probability_map(pm);

julia> q = fourier_map(pm);

julia> p_equivclasses = compute_equivalent_classes(p);

julia> p_equivclasses.parametrization
Dict{Tuple{Vararg{Int64}}, QQMPolyRingElem} with 5 entries:
  (1, 2, 1) => 1//4*a[1]*a[3]*b[2] + 1//4*a[2]*b[1]*b[3] + 1//2*b[1]*b[2]*b[3]
  (1, 1, 1) => 1//4*a[1]*a[2]*a[3] + 3//4*b[1]*b[2]*b[3]
  (1, 2, 2) => 1//4*a[1]*b[2]*b[3] + 1//4*a[2]*a[3]*b[1] + 1//2*b[1]*b[2]*b[3]
  (1, 2, 3) => 1//4*a[1]*b[2]*b[3] + 1//4*a[2]*b[1]*b[3] + 1//4*a[3]*b[1]*b[2] …
  (1, 1, 2) => 1//4*a[1]*a[2]*b[3] + 1//4*a[3]*b[1]*b[2] + 1//2*b[1]*b[2]*b[3]

julia> p_equivclasses.classes
Dict{Tuple{Vararg{Int64}}, Vector{Tuple{Vararg{Int64}}}} with 5 entries:
  (1, 2, 1) => [(1, 2, 1), (1, 3, 1), (1, 4, 1), (2, 1, 2), (2, 3, 2), (2, 4, 2…
  (1, 1, 1) => [(1, 1, 1), (2, 2, 2), (3, 3, 3), (4, 4, 4)]
  (1, 2, 2) => [(1, 2, 2), (1, 3, 3), (1, 4, 4), (2, 1, 1), (2, 3, 3), (2, 4, 4…
  (1, 2, 3) => [(1, 2, 3), (1, 2, 4), (1, 3, 2), (1, 3, 4), (1, 4, 2), (1, 4, 3…
  (1, 1, 2) => [(1, 1, 2), (1, 1, 3), (1, 1, 4), (2, 2, 1), (2, 2, 3), (2, 2, 4…

julia> q_equivclasses = compute_equivalent_classes(q);

julia> q_equivclasses.parametrization
Dict{Tuple{Vararg{Int64}}, QQMPolyRingElem} with 5 entries:
  (1, 1, 1) => x[1, 1]*x[2, 1]*x[3, 1]
  (2, 3, 4) => x[1, 2]*x[2, 2]*x[3, 2]
  (2, 2, 1) => x[3, 1]*x[1, 2]*x[2, 2]
  (1, 2, 2) => x[1, 1]*x[2, 2]*x[3, 2]
  (2, 1, 2) => x[2, 1]*x[1, 2]*x[3, 2]
```
"""
function compute_equivalent_classes(parametrization::Dict{Tuple{Vararg{Int64}}, QQMPolyRingElem})
  polys = unique(collect(values(parametrization)))
  polys = polys[findall(!is_zero, polys)]

  equivalent_keys = []
  for value in polys
      eqv_class = sort([key for key in keys(parametrization) if parametrization[key] == value])
      append!(equivalent_keys, [eqv_class])
  end

  param_classes = Dict{Tuple{Vararg{Int64}}, QQMPolyRingElem}(
                  k[1] => parametrization[k[1]] for k in equivalent_keys)

  classes = Dict{Tuple{Vararg{Int64}}, Vector{Tuple{Vararg{Int64}}}}(
            k[1] => k for k in equivalent_keys)
            
  return (parametrization=param_classes, classes=classes)
end

function equivalent_classes(pm::GroupBasedPhylogeneticModel; coordinates::Symbol=:fourier)
  if coordinates == :probabilities
    parametrization = probability_map(pm)
  elseif coordinates == :fourier
    parametrization = fourier_map(pm)
  end

  compute_equivalent_classes(parametrization).classes
end

function equivalent_classes(pm::PhylogeneticModel)
  parametrization = probability_map(pm)
  compute_equivalent_classes(parametrization).classes
end

@doc raw"""
    sum_equivalent_classes(equivalent_classes::NamedTuple{(:parametrization, :classes), Tuple{Dict{Tuple{Vararg{Int64}}, QQMPolyRingElem}, Dict{Tuple{Vararg{Int64}}, Vector{Tuple{Vararg{Int64}}}}}})

Take the output of the function `compute_equivalent_classes` for `PhylogeneticModel` and multiply by a factor to obtain probabilities as specified on the original small trees database.

# Examples
```jldoctest
julia> pm = jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]));

julia> p = probability_map(pm);

julia> p_equivclasses = compute_equivalent_classes(p);

julia> sum_equivalent_classes(p_equivclasses)
Dict{Tuple{Int64, Int64, Int64}, QQMPolyRingElem} with 5 entries:
  (1, 2, 1) => 3*a[1]*a[3]*b[2] + 3*a[2]*b[1]*b[3] + 6*b[1]*b[2]*b[3]
  (1, 1, 1) => a[1]*a[2]*a[3] + 3*b[1]*b[2]*b[3]
  (1, 2, 2) => 3*a[1]*b[2]*b[3] + 3*a[2]*a[3]*b[1] + 6*b[1]*b[2]*b[3]
  (1, 2, 3) => 6*a[1]*b[2]*b[3] + 6*a[2]*b[1]*b[3] + 6*a[3]*b[1]*b[2] + 6*b[1]*…
  (1, 1, 2) => 3*a[1]*a[2]*b[3] + 3*a[3]*b[1]*b[2] + 6*b[1]*b[2]*b[3]
```
"""
function sum_equivalent_classes(equivalent_classes::NamedTuple{(:parametrization, :classes), Tuple{Dict{Tuple{Vararg{Int64}}, QQMPolyRingElem}, Dict{Tuple{Vararg{Int64}}, Vector{Tuple{Vararg{Int64}}}}}})
  param = equivalent_classes.parametrization
  classes = equivalent_classes.classes
  return Dict(key => param[key]*length(classes[key]) for key in keys(param))
end


##############################################
#### SPECIALIZED FOURIER TRANSFORM MATRIX ####
##############################################

@doc raw"""
    specialized_fourier_transform(pm::GroupBasedPhylogeneticModel, p_classes::Dict{Tuple{Vararg{Int64}}, Vector{Tuple{Vararg{Int64}}}}, q_classes::Dict{Tuple{Vararg{Int64}}, Vector{Tuple{Vararg{Int64}}}})

Reparametrize between a model specification in terms of probability and Fourier coordinates. The input of equivalent classes is optional, if they are not entered they will be computed.

# Examples
```jldoctest
julia> pm = jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]));

julia> p_equivclasses = compute_equivalent_classes(probability_map(pm));

julia> q_equivclasses = compute_equivalent_classes(fourier_map(pm));

julia> specialized_fourier_transform(pm, p_equivclasses.classes, q_equivclasses.classes)
5×5 Matrix{QQMPolyRingElem}:
 1  1      1      1      1
 1  -1//3  -1//3  1      -1//3
 1  -1//3  1      -1//3  -1//3
 1  1      -1//3  -1//3  -1//3
 1  -1//3  -1//3  -1//3  1//3
```
"""
function specialized_fourier_transform(pm::GroupBasedPhylogeneticModel, p_classes::Dict{Tuple{Vararg{Int64}}, Vector{Tuple{Vararg{Int64}}}}, 
                                       q_classes::Dict{Tuple{Vararg{Int64}}, Vector{Tuple{Vararg{Int64}}}})
    R = probability_ring(pm)
    ns = number_states(pm)
    
    np = length(p_classes)
    nq = length(q_classes)
    
    ## We need to sort the equivalence classes: both inside each class as well as the collection of classes. 
    keys_p_classes = collect(keys(p_classes))
    sort!(keys_p_classes)
    keys_q_classes = collect(keys(q_classes))
    sort!(keys_q_classes)
    
    H = R.(hadamard(matrix_space(ZZ, ns, ns)))
    
    specialized_ft_matrix = R.(Int.(zeros(nq, np)))
    for i in 1:nq
      current_fourier_class = q_classes[keys_q_classes[i]]
      for j in 1:np
        current_prob_class = p_classes[keys_p_classes[j]]
        current_entriesin_M = [prod([H[y,x] for (x,y) in zip(p,q)]) for p in current_prob_class, q in current_fourier_class]
        specialized_ft_matrix[i,j] = R.(1//(length(current_prob_class)*length(current_fourier_class))*sum(current_entriesin_M))
      end
    end
    return specialized_ft_matrix
end

function specialized_fourier_transform(pm::GroupBasedPhylogeneticModel)
  p_equivclasses = compute_equivalent_classes(probability_map(pm))
  q_equivclasses = compute_equivalent_classes(fourier_map(pm))
  specialized_fourier_transform(pm, p_equivclasses.classes, q_equivclasses.classes)
end

@doc raw"""
    inverse_specialized_fourier_transform(pm::GroupBasedPhylogeneticModel, p_classes::Dict{Tuple{Vararg{Int64}}, Vector{Tuple{Vararg{Int64}}}}, q_classes::Dict{Tuple{Vararg{Int64}}, Vector{Tuple{Vararg{Int64}}}})

Reparametrize between a model specification in terms of Fourier and probability coordinates.

# Examples
```jldoctest
julia> pm = jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]));

julia> p_equivclasses = compute_equivalent_classes(probability_map(pm));

julia> q_equivclasses = compute_equivalent_classes(fourier_map(pm));

julia> inverse_specialized_fourier_transform(pm, p_equivclasses.classes, q_equivclasses.classes)
5×5 Matrix{QQMPolyRingElem}:
 1//16  3//16   3//16   3//16   3//8
 3//16  -3//16  -3//16  9//16   -3//8
 3//16  -3//16  9//16   -3//16  -3//8
 3//16  9//16   -3//16  -3//16  -3//8
 3//8   -3//8   -3//8   -3//8   3//4
```
"""
function inverse_specialized_fourier_transform(pm::GroupBasedPhylogeneticModel, p_classes::Dict{Tuple{Vararg{Int64}}, Vector{Tuple{Vararg{Int64}}}}, 
                                               q_classes::Dict{Tuple{Vararg{Int64}}, Vector{Tuple{Vararg{Int64}}}})
  R = probability_ring(pm)
  ns = number_states(pm)

  np = length(p_classes)
  nq = length(q_classes)

  ## We need to sort the equivalence classes: both inside each class as well as the collection of classes. 
  keys_p_classes = collect(keys(p_classes))
  sort!(keys_p_classes)
  keys_q_classes = collect(keys(q_classes))
  sort!(keys_q_classes)

  H = R.(hadamard(matrix_space(ZZ, ns, ns)))
  Hinv = 1//ns * H 

  inverse_spec_ft_matrix = R.(Int.(zeros(np, nq)))
  for i in 1:np
    current_prob_class = p_classes[keys_p_classes[i]]
    for j in 1:nq
      current_fourier_class = q_classes[keys_q_classes[j]]
      current_entriesin_Minv = [prod([Hinv[x,y] for (x,y) in zip(p,q)]) for p in current_prob_class, q in current_fourier_class] 
      inverse_spec_ft_matrix[i,j] = R.(sum(current_entriesin_Minv))
    end
  end
  return inverse_spec_ft_matrix
end

function inverse_specialized_fourier_transform(pm::GroupBasedPhylogeneticModel)
  p_equivclasses = compute_equivalent_classes(probability_map(pm))
  q_equivclasses = compute_equivalent_classes(fourier_map(pm))
  inverse_specialized_fourier_transform(pm, p_equivclasses.classes, q_equivclasses.classes)
end
