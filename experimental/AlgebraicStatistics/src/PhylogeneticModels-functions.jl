
### GRAPH FUNCTIONS ###
# ---------------------

function interior_nodes(graph::Graph)
  big_graph = Polymake.graph.Graph(ADJACENCY = pm_object(graph))
  degrees = big_graph.NODE_DEGREES
  return findall(>(1), degrees)
end

function leaves(graph::Graph)
  big_graph = Polymake.graph.Graph(ADJACENCY = pm_object(graph))
  degrees = big_graph.NODE_DEGREES
  return findall(==(1), degrees)
end

function n_leaves(graph::Graph)
  length(leaves(graph))
end

function root(graph::Graph)
  n_parents = [length(inneighbors(graph, v)) for v in 1:n_vertices(graph)]
  return findall(x -> x == 0, n_parents)[1]
end

function sort_edges(graph::Graph)
  edgs = collect(edges(graph))
  leaves_idx = findall(edge -> dst(edge) in Oscar.leaves(graph), edgs)
  return edgs[vcat(leaves_idx, setdiff(1:length(edgs), leaves_idx))]
end


### PARAMETRIZATIONS ###
# ----------------------

function leaves_indices(PM::PhylogeneticModel)
  leave_nodes = leaves(graph(PM))
  leaves_indices = collect(Iterators.product(
    [tuple(1:n_states(PM)...) for _ in leave_nodes]...))

  return leaves_indices
end

function fully_observed_probability(PM::PhylogeneticModel, vertices_states::Dict{Int, Int})
  gr = graph(PM)
  root_dist = root_distribution(PM)

  r = root(gr)
  monomial = root_dist[vertices_states[r]]
  
  for edge in edges(gr)
    state_parent = vertices_states[src(edge)]
    state_child = vertices_states[dst(edge)]
    # get the symbolfrom the transition matrix signature
    sym = entry_transition_matrix(PM, edge, state_parent, state_child)
    # we can try and avoid this here if this becomes a bottle neck
    # it's related to the comment below about using a polynomial context
    # i.e., this is just the adding of exponents which are integer vectors
    # which would be much faster than polynomial multiplication
    monomial = monomial * sym
  end

  return monomial
end

function leaves_probability(PM::PhylogeneticModel, leaves_states::Dict{Int, Int})
  gr = graph(PM)
  int_nodes = interior_nodes(gr)

  interior_indices = collect.(Iterators.product([collect(1:n_states(PM)) for _ in int_nodes]...))  
  vertices_states = leaves_states

  poly = 0
  # Might be useful in the future to use a polynomial ring context
  for labels in interior_indices
    for (int_node, label) in zip(int_nodes, labels)
      vertices_states[int_node] = label
    end
    poly = poly + fully_observed_probability(PM, vertices_states)
  end 
  return poly
end 



### OTHER AUXILIARY FUNCTIONS ###
# -------------------------------

# I would like to get rid of this! Is there a better way?
function index(x::QQMPolyRingElem)

    m = match(r"\[(.*)\]", string(x))
    idx = split(m.captures[1], ',')

    [parse(Int, n) for n in idx if !isempty(n)]
end

