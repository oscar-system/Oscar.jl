# Add your new types, functions, and methods here.

struct PhylogeneticModel
  n_states::Int
  graph::Graph{Directed}
  prob_ring::MPolyRing{QQFieldElem}
  fourier_ring::MPolyRing{QQFieldElem}
  trans_matrices::Dict{Edge, MatElem{QQMPolyRingElem}}
  fourier_params::Dict{Edge, Vector{QQMPolyRingElem}}
  group::Vector{Vector{Int64}}
end

graph(pm::PhylogeneticModel) = pm.graph
transition_matrices(pm::PhylogeneticModel) = pm.trans_matrices
number_states(pm::PhylogeneticModel) = pm.n_states
probablities_ring(pm::PhylogeneticModel) = pm.prob_ring
fourier_ring(pm::PhylogeneticModel) = pm.fourier_ring
fourier_parameters(pm::PhylogeneticModel) = pm.fourier_params
group_model(pm::PhylogeneticModel) = pm.group

function jukes_cantor_model(graph::Graph{Directed})
  ns = 4
  ne = n_edges(graph)
  R, list_a, list_b = polynomial_ring(QQ, :a => 1:ne, :b => 1:ne)
  
  matrices = Dict{Edge, MatElem}(e => matrix(R, [
    a b b b
    b a b b
    b b a b
    b b b a]) for (a,b,e) in zip(list_a, list_b, edges(graph))
  )

  S, list_x = polynomial_ring(QQ, :x => (1:ne, 1:2))
  fourier_param = Dict{Edge, Vector{QQMPolyRingElem}}(e => 
          [list_x[i,1], list_x[i,2], list_x[i,2], list_x[i,2]] for (i, e) in zip(1:ne, edges(graph)))
  
  group = [[0,0], [0,1], [1,0], [1,1]]

  return PhylogeneticModel(ns, graph, R, S, matrices, fourier_param, group)
end

function interior_nodes(graph::Graph)
  big_graph = Polymake.graph.Graph(ADJACENCY = pm_object(graph))
  degrees = big_graph.NODE_DEGREES
  return findall(x -> x > 1, degrees)
end

function leaves(graph::Graph)
  big_graph = Polymake.graph.Graph(ADJACENCY = pm_object(graph))
  degrees = big_graph.NODE_DEGREES
  return findall(x -> x == 1, degrees)
end


function vertex_descendants(v::Int, gr::Graph, desc::Vector{Any})
  lvs = leaves(gr)
  outn = outneighbors(gr, v)

  if v in lvs
    return([v])
  end

  innodes = setdiff(outn, lvs)
  d = unique(append!(desc, intersect(outn, lvs)))
 
  if length(innodes) > 0
      for i in innodes
          d = vertex_descendants(i, gr, d)
      end
      return(d)
  end

  return(d)
end



function monomial_parametrization(pm::PhylogeneticModel, states::Dict{Int, Int})
  gr = graph(pm)
  tr_mat = transition_matrices(pm)

  monomial = 1
  for edge in edges(gr)
      stateParent = states[src(edge)]
      stateChild = states[dst(edge)]
      monomial = monomial * tr_mat[edge][stateParent, stateChild]
  end

  return(monomial)
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

function probability_map(pm::PhylogeneticModel)
  lvs_nodes = leaves(graph(pm))
  n_states = number_states(pm)

  leaves_indices = collect.(Iterators.product([collect(1:n_states) for _ in lvs_nodes]...))
  probability_coordinates = Dict(leaves_states => probability_parametrization(pm, leaves_states) for leaves_states in leaves_indices)
  return(probability_coordinates)
end

#### SPECIALIZED FOURIER TRANSFORM MATRIX ####

function hadamardmatrix()
  H = [1 1 1 1
       1 1 -1 -1 
       1 -1 1 -1
       1 -1 -1 1]
  return H
end

# We need the parameter leaf_number in order to determine the right size of
# the Kronecker product. 
# TODO: line 138 mutates f_equivclasses, i.e. an object from outside the function.
# This makes the function a !-function. We do not want f_equivclasses to be changed, 
# so I add the missing equivalence class again in line 175. 
# We should find out whether there is a better way to do this.
function specialized_fourier_transform(pm::PhylogeneticModel, p_equivclasses::Dict{Vector{Vector{Int64}}, QQMPolyRingElem},f_equivclasses::Dict{Vector{Vector{Int64}}, QQMPolyRingElem})
  R = polyn_ring(pm)

  class0 = findall(x -> x ==0, f_equivclasses)[1]
  delete!(f_equivclasses, class0)

  np = length(p_equivclasses)
  nq = length(f_equivclasses)

  ## We need to sort the equivalence classes: both inside each class as well as the collection of classes. 
  p_equivclasses_sorted = collect(keys(p_equivclasses))
  for p_eqclass in p_equivclasses_sorted
      sort!(p_eqclass)
  end
  sort!(p_equivclasses_sorted)

  f_equivclasses_sorted = collect(keys(f_equivclasses))
  for f_eqclass in f_equivclasses_sorted
      sort!(f_eqclass)
  end
  sort!(f_equivclasses_sorted)

  H = hadamardmatrix()

  specialized_ft_matrix = R.(Int.(zeros(nq, np)))
  for i in 1:nq
      current_fourier_classes = f_equivclasses_sorted[i]
      for j in 1:np
          current_prob_classes = p_equivclasses_sorted[j]
          current_entriesin_M = [prod([H[y,x] for (x,y) in zip(p,q)]) for p in current_prob_classes, q in current_fourier_classes]
          specialized_ft_matrix[i,j] = R.(1//(length(current_prob_classes)*length(current_fourier_classes))*sum(current_entriesin_M))
      end
  end
  get!(f_equivclasses, class0, R(0))
  return specialized_ft_matrix
end

function inverse_specialized_fourier_transform(pm::PhylogeneticModel, p_equivclasses::Dict{Vector{Vector{Int64}}, QQMPolyRingElem},f_equivclasses::Dict{Vector{Vector{Int64}}, QQMPolyRingElem})
  R = polyn_ring(pm)

  class0 = findall(x -> x ==0, f_equivclasses)[1]
  delete!(f_equivclasses, class0)

  np = length(p_equivclasses)
  nq = length(f_equivclasses)

  ## We need to sort the equivalence classes: both inside each class as well as the collection of classes. 
  p_equivclasses_sorted = collect(keys(p_equivclasses))
  for p_eqclass in p_equivclasses_sorted
      sort!(p_eqclass)
  end
  sort!(p_equivclasses_sorted)

  f_equivclasses_sorted = collect(keys(f_equivclasses))
  for f_eqclass in f_equivclasses_sorted
      sort!(f_eqclass)
  end
  sort!(f_equivclasses_sorted)

  H = hadamardmatrix()
  Hinv = 1//4 * H 

  inverse_spec_ft_matrix = R.(Int.(zeros(np, nq)))
  for i in 1:np
      current_prob_class = p_equivclasses_sorted[i]
      for j in 1:nq
          current_fourier_class = f_equivclasses_sorted[j]
          current_entriesin_Minv = [prod([Hinv[x,y] for (x,y) in zip(p,q)]) for p in current_prob_class, q in current_fourier_class] 
          inverse_spec_ft_matrix[i,j] = R.(sum(current_entriesin_Minv))
      end
  end
  get!(f_equivclasses, class0, R(0))
  return inverse_spec_ft_matrix
end

-
# @doc raw"""
#     my_access_func(S::ExampleStruct)

# This is a dummy sample function to teach the use of docstrings.
# """

function group_sum(pm::PhylogeneticModel, states::Vector{Int})
  group = group_model(pm)

  if (length(states) == 1)
    return(group[states[1]])
  end

  sum(group[states]).%2
end


function which_group_element(pm::PhylogeneticModel, elem::Vector{Vector{Int64}})
  group = group_model(pm)
  return findall([all(group[i].==elem) for i in 1:length(group)])[1]
end

function monomial_fourier(pm::PhylogeneticModel, leaves_states::Vector{Int})
  gr = graph(pm)
  param = fourier_parameters(pm)
  monomial = 1
  for edge in edges(gr)
    dsc = vertex_descendants(dst(edge), gr, [])
    elem = group_sum(pm, leaves_states[dsc])
    monomial = monomial * param[edge][which_group_element(pm, elem)]
  end

  return(monomial)
end

function fourier_parametrization(pm::PhylogeneticModel, leaves_states::Vector{Int})
  S = fourier_ring(pm)
  if group_sum(pm, leaves_states) == [0,0]
    poly = monomial_fourier(pm, leaves_states)
  else 
    poly = S(0)
  end

  return poly
end 

function fourier_map(pm::PhylogeneticModel)
  lvs_nodes = leaves(graph(pm))
  n_states = number_states(pm)

  leaves_indices = collect.(Iterators.product([collect(1:n_states) for _ in lvs_nodes]...))
  fourier_coordinates = Dict(leaves_states => fourier_parametrization(pm, leaves_states) for leaves_states in leaves_indices)
  return(fourier_coordinates)
end
