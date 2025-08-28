##########################
#### GROUP OPERATIONS ####
##########################

function group_sum(pm::GroupBasedPhylogeneticModel, states::Vector{Int})
  group = group_of_model(pm)
  return sum(group[states])
end

function is_zero_group_sum(pm::GroupBasedPhylogeneticModel, states::Vector{Int})
  return group_sum(pm, states) == zero(parent(group_of_model(pm)[1]))
end
  
function which_group_element(pm::GroupBasedPhylogeneticModel, elem::FinGenAbGroupElem)
  group = group_of_model(pm)
  return findall([all(g==elem) for g in group])[1]
end


########################################
#### AUXILIARY FUNCTIONS FOR GRAPHS ####
########################################

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

function root(graph::Graph)
  n_parents = [length(inneighbors(graph, v)) for v in 1:n_vertices(graph)]
  return findfirst(x -> x == 0, n_parents)[1]
end


########################
#### SORT FUNCTIONS ####
########################

function sort_edges(graph::Graph)
  edgs = collect(edges(graph))
  leaves_idx = findall(edge -> dst(edge) in Oscar.leaves(graph), edgs)
  return edgs[vcat(leaves_idx, setdiff(1:length(edgs), leaves_idx))]
end
