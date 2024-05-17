##########################
#### GROUP OPERATIONS ####
##########################
# This function implements the group operation on Z2 x Z2 and is subsequently 
# used to check the condition given in the computation of the Fourier
# coordinates and to calculate the sum of the leaves beneath a given edge. 

function group_sum(pm::GroupBasedPhylogeneticModel, states::Vector{Int})
  group = group_of_model(pm)
  return sum(group[states]).%2
end

function is_zero_group_sum(pm::GroupBasedPhylogeneticModel, states::Vector{Int})
  ng = length(states)
  return group_sum(pm, [states[1]]) == group_sum(pm, states[2:ng])
end
  
function which_group_element(pm::GroupBasedPhylogeneticModel, elem::Vector{Int64})
  group = group_of_model(pm)
  return findall([all(group[i].==elem) for i in 1:length(group)])[1]
end


########################################
#### AUXILIARY FUNCTIONS FOR GRAPHS ####
########################################

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
    return [v]
  end

  innodes = setdiff(outn, lvs)
  d = unique(append!(desc, intersect(outn, lvs)))
  
  if length(innodes) > 0
    for i in innodes
      d = vertex_descendants(i, gr, d)
    end
    return d 
  end

  return d
end

function cherries(graph::Graph)
  lvs = leaves(graph)
  cherry = []
  for l in lvs
    in_node = inneighbors(graph,l)[1]
    lvs_cherr = outneighbors(graph, inneighbors(graph,l)[1])
    if issubset(lvs_cherr, lvs) == 2
      cherry = append!(cherry, [[Edge(in_node, lvs_cherr[1]), Edge(in_node, lvs_cherr[2])]])
    end
  end
  
  return unique(cherry)
end

function root(graph::Graph)
  n_parents = [length(inneighbors(graph, v)) for v in 1:n_vertices(graph)]
  return findall(x -> x == 0, n_parents)[1]
end


############################
#### ORDERING FUNCTIONS ####
############################

function order_edges(graph::Graph)
  edgs = collect(edges(graph))
  leaves_idx = findall(edge -> dst(edge) in Oscar.leaves(graph), edgs)
  return edgs[vcat(leaves_idx, setdiff(1:length(edgs), leaves_idx))]
end
