###################################################################################
#
#       Auxiliary graph functions
#
###################################################################################

interior_nodes(graph::Graph) = findall(>(1), outdegree(graph))

leaves(graph::Graph) = findall(iszero, outdegree(graph))
n_leaves(graph::Graph) = length(leaves(graph))
n_leaves(pt::PhylogeneticTree) = n_leaves(adjacency_tree(pt))


function root(graph::Graph)
  n_parents = [length(inneighbors(graph, v)) for v in 1:n_vertices(graph)]
  return findall(iszero, n_parents)[1]
end

function sort_edges(graph::Graph)
  edgs = collect(edges(graph))
  leaves_idx = findall(edge -> dst(edge) in Oscar.leaves(graph), edgs)
  return edgs[vcat(leaves_idx, setdiff(1:length(edgs), leaves_idx))]
end

sort_edges(pt::PhylogeneticTree) = sort_edges(adjacency_tree(pt))

function vertex_descendants(gr::Graph{Directed}, v::Int, desc::Vector{Any} = [])
  lvs = leaves(gr)
  outn = outneighbors(gr, v)
  if v in lvs
    return [v]
  end

  innodes = setdiff(outn, lvs)
  d = unique(append!(desc, intersect(outn, lvs)))
  
  if length(innodes) > 0
    for i in innodes
      d = vertex_descendants(gr, i, d)
    end
    return d 
  end

  return d
end

vertex_descendants(pt::PhylogeneticTree, v::Int, desc::Vector{Any} = []) = vertex_descendants(adjacency_tree(pt), v, desc)


###################################################################################
#
#       Auxiliary functions to access parameters
#
###################################################################################


function entry_transition_matrix(PM::PhylogeneticModel{<:AbstractGraph{Directed}, <: VarName, L, T}, i::Int, j::Int, e::Edge) where {L, T}
  tr_mat = transition_matrix(PM)
  parameter_ring(PM)[2][tr_mat[i,j], e]
end

function entry_transition_matrix(PM::PhylogeneticModel{<:AbstractGraph{Directed}, <: VarName, L, T}, i::Int, j::Int, u::Int, v::Int) where {L, T}
  tr_mat = transition_matrix(PM)
  parameter_ring(PM)[2][tr_mat[i,j], Edge(u,v)]
end

function entry_transition_matrix(PM::PhylogeneticModel{<:AbstractGraph{Directed}, <: MPolyRingElem, L, <: T}, i::Int, j::Int, e::Edge) where {L, T}
  tr_mat = transition_matrix(PM)
  parameter_ring(PM)[2][e](tr_mat[i,j])
end

function entry_transition_matrix(PM::PhylogeneticModel{<:AbstractGraph{Directed}, <: MPolyRingElem, L, T}, i::Int, j::Int, u::Int, v::Int) where {L, T}
  tr_mat = transition_matrix(PM)
  parameter_ring(PM)[2][Edge(u,v)](tr_mat[i,j])
end


# Is this fine or entry_transition_matrix should only be defined for a PhyloModel?
function entry_transition_matrix(PM::GroupBasedPhylogeneticModel, i::Int, j::Int, e::Edge)
  entry_transition_matrix(phylogenetic_model(PM), i, j, e)
end

# ?
function entry_transition_matrix(PM::GroupBasedPhylogeneticModel, i::Int, j::Int, u::Int, v::Int)
  entry_transition_matrix(phylogenetic_model(PM), i, j, u, v)
end

function entry_root_distribution(PM::PhylogeneticModel, i::Int)
  parameter_ring(PM)[3][i]
end

#?
function entry_root_distribution(PM::GroupBasedPhylogeneticModel, i::Int)
  entry_root_distribution(phylogenetic_model(PM), i)
  
end

function entry_fourier_parameter(PM::GroupBasedPhylogeneticModel, i::Int, e::Edge)
  x = fourier_parameters(PM)
  parameter_ring(PM)[2][x[i], e]
end

function entry_fourier_parameter(PM::GroupBasedPhylogeneticModel, i::Int, u::Int, v::Int)
  x = fourier_parameters(PM)
  parameter_ring(PM)[2][x[i], Edge(u,v)]
end


###################################################################################
#
#       Auxiliary functions to compute the parametrizations
#
###################################################################################


function leaves_indices(PM::PhylogeneticModel)
  leave_nodes = leaves(graph(PM))
  leaves_indices = collect(Iterators.product(
    [tuple(1:n_states(PM)...) for _ in leave_nodes]...))

  return leaves_indices
end
function leaves_indices(PM::GroupBasedPhylogeneticModel)
  return leaves_indices(phylogenetic_model(PM))
end

function fully_observed_probability(PM::PhylogeneticModel, vertices_states::Dict{Int, Int})
  gr = graph(PM)

  r = root(gr)
  monomial = entry_root_distribution(PM, vertices_states[r])
  
  for edge in edges(gr)
    state_parent = vertices_states[src(edge)]
    state_child = vertices_states[dst(edge)]
    # get the symbolfrom the transition matrix signature
    sym = entry_transition_matrix(PM, state_parent, state_child, edge)
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

function monomial_fourier(PM::GroupBasedPhylogeneticModel, leaves_states::Dict{Int, Int})
  gr = graph(PM)
  monomial = 1
  
  for edge in edges(gr)
    dsc = vertex_descendants(gr, dst(edge))
    dsc_states = filter(pair -> pair.first in dsc, leaves_states)
    group_elem = group_sum(PM, dsc_states)
    monomial = monomial * entry_fourier_parameter(PM, which_group_element(PM, group_elem), edge)
  end
  return monomial
end

function leaves_fourier(PM::GroupBasedPhylogeneticModel, leaves_states::Dict{Int, Int})
  
  if is_zero_group_sum(PM, leaves_states) 
    return monomial_fourier(PM, leaves_states)
  end 

  S = parameter_ring(PM)[1]
  return S(0)
  
end 
  

###################################################################################
#
#       Auxiliary group operation functions
#
###################################################################################


function group_sum(PM::GroupBasedPhylogeneticModel, states::Dict{Int, Int})
  G = group(PM)
  return sum(G[collect(values(states))])
end

is_zero_group_sum(PM::GroupBasedPhylogeneticModel, states::Dict{Int, Int}) =  iszero(group_sum(PM, states))


function which_group_element(PM::GroupBasedPhylogeneticModel, elem::FinGenAbGroupElem)
  G = group(PM)
  return findall([all(g==elem) for g in G])[1]
end
