# Add your new types, functions, and methods here.
# :)

struct PhylogeneticModel
  n_states::Int
  graph::Graph{Directed}
  R::MPolyRing{QQFieldElem}
  trans_matrices::Dict{Edge, MatElem{QQMPolyRingElem}}
end

graph(pm::PhylogeneticModel) = pm.graph
transition_matrices(pm::PhylogeneticModel) = pm.trans_matrices
number_states(pm::PhylogeneticModel) = pm.n_states

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
  return PhylogeneticModel(ns, graph, R, matrices)
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

end

-
# @doc raw"""
#     my_access_func(S::ExampleStruct)

# This is a dummy sample function to teach the use of docstrings.
# """
