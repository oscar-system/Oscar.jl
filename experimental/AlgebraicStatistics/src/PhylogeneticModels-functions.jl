###################################################################################
#
#       Auxiliary graph functions
#
###################################################################################

interior_nodes(graph::Graph) = findall(>=(1), outdegree(graph))

leaves(graph::Graph) = findall(iszero, outdegree(graph))
n_leaves(graph::Graph) = length(leaves(graph))
n_leaves(pt::PhylogeneticTree) = n_leaves(adjacency_tree(pt))


function roots(graph::Graph)
  n_parents = [length(inneighbors(graph, v)) for v in 1:n_vertices(graph)]
  return findall(iszero, n_parents)
end

root(graph::Graph) = roots(graph)[1]


function sort_edges(graph::Graph)
  edgs = collect(edges(graph))
  leaves_idx = findall(edge -> dst(edge) in leaves(graph), edgs)
  return edgs[vcat(leaves_idx, setdiff(1:length(edgs), leaves_idx))]
end

sort_edges(N::PhylogeneticNetwork) = sort_edges(graph(N))
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

function biconnected_components(g::Graph{Undirected})
    im = Polymake.call_function(:graph, :biconnected_components, pm_object(g))::IncidenceMatrix
    return [Vector(Polymake.row(im,i)) for i in 1:Polymake.nrows(im)]
end

###################################################################################
#
#       Auxiliary graph functions for phylogenetic networks
#
###################################################################################

function is_phylogenetic_network(G::Graph{Directed})
    ## We assume phylogenetic networks are binary

    A = matrix(ZZ, adjacency_matrix(G))
    indegrees = sum(A[i,:] for i in 1:size(A)[1])
    outdegrees = sum(A[:,i] for i in 1:size(A)[1])

    roots = findall(indegrees .== 0) #findall((indegrees .== 0) .& (outdegrees .== 2))
    lvs = findall((indegrees .== 1) .& (outdegrees .== 0))
    tree_nodes = findall((indegrees .== 1) .& (outdegrees .<= 2) .& (outdegrees .> 0))
    h_nodes = findall((indegrees .== 2) .& (outdegrees .== 1))

    if length(roots) != 1
        return false
    end

    if length(lvs) == 0
        return false
    end

    if length(vertices(G)) != length(vcat(roots, lvs, tree_nodes, h_nodes))
        return false
    end

    return true
end

function hybrid_vertices(G::Graph{Directed})
    A = matrix(ZZ, adjacency_matrix(G))
    indegrees = sum(A[i,:] for i in 1:size(A)[1])
    outdegrees = sum(A[:,i] for i in 1:size(A)[1])

    findall((indegrees .== 2) .& (outdegrees .== 1))
end

function hybrid_edges(G::Graph{Directed})
    hybrid_nodes = hybrid_vertices(G)
    [[Edge(j, i) for j in sort(inneighbors(G, i))] for i in hybrid_nodes]
end

function hybrid_edges(G::Graph{Directed}, i)
    hybrid_nodes = hybrid_vertices(G)
    if !(i in hybrid_nodes)
        error("$i is not a hybrid node.")
    end
    [Edge(j, i) for j in sort(inneighbors(G, i))]
end

function level_phylogenetic_network(G::Graph{Directed})
    if !is_phylogenetic_network(G)
        error("$G is not a phylogenetic network.")
    end

    h_nodes = hybrid_vertices(G)
    bicon_comps = biconnected_components(graph_from_edges(Undirected,edges(G)))

    maximum([length(intersect(component, h_nodes)) for component in bicon_comps])
end

function tree_edges(N::PhylogeneticNetwork)
  hyb = hybrids(N)
  egdes_N = collect(edges(N))
  egdes_N[findall(e -> !(dst(e) in collect(keys(hyb))), egdes_N)]
end

###################################################################################
#
#       Auxiliary functions to access parameters
#
###################################################################################


function entry_transition_matrix(PM::PhylogeneticModel{<:AbstractGraph{Directed}, L, <: VarName, T}, i::Int, j::Int, e::Edge) where {L, T}
  tr_mat = transition_matrix(PM)
  parameter_ring(PM)[2][tr_mat[i,j], e]
end

function entry_transition_matrix(PM::PhylogeneticModel{<:AbstractGraph{Directed}, L, <: VarName, T}, i::Int, j::Int, u::Int, v::Int) where {L, T}
  entry_transition_matrix(PM, i, j, Edge(u,v))
end

function entry_transition_matrix(PM::PhylogeneticModel{<:AbstractGraph{Directed}, L, <: MPolyRingElem, <: T}, i::Int, j::Int, e::Edge) where {L, T}
  tr_mat = transition_matrix(PM)
  parameter_ring(PM)[2][e](tr_mat[i,j])
end

function entry_transition_matrix(PM::PhylogeneticModel{<:AbstractGraph{Directed}, L, <: MPolyRingElem, T}, i::Int, j::Int, u::Int, v::Int) where {L, T}
  entry_transition_matrix(PM, i, j, Edge(u,v))
end

# Is this fine or entry_transition_matrix should only be defined for a PhyloModel?
entry_transition_matrix(PM::GroupBasedPhylogeneticModel, i::Int, j::Int, e::Edge) = entry_transition_matrix(phylogenetic_model(PM), i, j, e)
entry_transition_matrix(PM::GroupBasedPhylogeneticModel, i::Int, j::Int, u::Int, v::Int) = entry_transition_matrix(phylogenetic_model(PM), i, j, u, v)

function entry_root_distribution(PM::PhylogeneticModel, i::Int)
  parameter_ring(PM)[3][i]
end

entry_root_distribution(PM::GroupBasedPhylogeneticModel, i::Int) = entry_root_distribution(phylogenetic_model(PM), i)

function entry_fourier_parameter(PM::GroupBasedPhylogeneticModel, i::Int, e::Edge)
  x = fourier_parameters(PM)
  parameter_ring(PM)[2][x[i], e]
end

entry_fourier_parameter(PM::GroupBasedPhylogeneticModel, i::Int, u::Int, v::Int) = entry_fourier_parameter(PM, i, Edge(u,v))

function entry_hybrid_parameter(PM::PhylogeneticModel{<: PhylogeneticNetwork}, e::Edge)
  if !(dst(e) in collect(keys(hybrids(graph(PM)))))
    error("Edge $e id not a hybrid edge in the Phylogenetic network $(graph(PM))")
  end

  parameter_ring(PM)[4][e]
end

function entry_hybrid_parameter(PM::GroupBasedPhylogeneticModel{<: PhylogeneticNetwork}, e::Edge) 
  if !(dst(e) in collect(keys(hybrids(graph(PM)))))
    error("Edge $e id not a hybrid edge in the Phylogenetic network $(graph(PM))")
  end

  parameter_ring(PM)[3][e]
end

function entry_hybrid_parameter(PM::Union{GroupBasedPhylogeneticModel, PhylogeneticModel}, u::Int, v::Int)
  entry_hybrid_parameter(PM, Edge(u,v))
end

###################################################################################
#
#       Auxiliary functions to compute the parametrizations
#
###################################################################################

root(PM::PhylogeneticModel) = _root(adjacency_tree(graph(PM)))

function leaves_indices(PM::PhylogeneticModel)
  leave_nodes = leaves(graph(PM))
  leaves_indices = collect(Iterators.product(
    [tuple(1:n_states(PM)...) for _ in leave_nodes]...))

  return leaves_indices
end

leaves_indices(PM::GroupBasedPhylogeneticModel) = leaves_indices(phylogenetic_model(PM))

function hybrid_indices(PM::Union{GroupBasedPhylogeneticModel{GT}, PhylogeneticModel{GT}}) where {GT <: PhylogeneticNetwork} 
  hyb = hybrids(graph(PM))
  hyb_indices = collect(Iterators.product([tuple(1:2...) for _ in keys(hyb)]...))

  return hyb_indices
end

function fully_observed_probability(PM::PhylogeneticModel{<:PhylogeneticTree}, vertices_states::Dict{Int, Int}) 
  gr = graph(PM)

  r = root(PM)
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

function fully_observed_probability(PM::PhylogeneticModel{<:PhylogeneticNetwork}, vertices_states::Dict{Int, Int}, subtree::Graph{Directed}) 

  r = root(subtree)
  monomial = entry_root_distribution(PM, vertices_states[r])
  
  for edge in edges(subtree)
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

function leaves_probability(PM::PhylogeneticModel{<:PhylogeneticTree}, leaves_states::Dict{Int, Int})
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

function leaves_probability(PM::PhylogeneticModel{<:PhylogeneticNetwork}, leaves_states::Dict{Int, Int}, subtree::Graph{Directed})
  int_nodes = interior_nodes(subtree)

  interior_indices = collect.(Iterators.product([collect(1:n_states(PM)) for _ in int_nodes]...))  
  vertices_states = leaves_states

  poly = 0
  # Might be useful in the future to use a polynomial ring context
  for labels in interior_indices
    for (int_node, label) in zip(int_nodes, labels)
      vertices_states[int_node] = label
    end
    poly = poly + fully_observed_probability(PM, vertices_states, subtree)
  end 
  return poly
end 

function monomial_fourier(PM::GroupBasedPhylogeneticModel{<:PhylogeneticTree}, leaves_states::Dict{Int, Int})
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

function monomial_fourier(PM::GroupBasedPhylogeneticModel{<:PhylogeneticNetwork}, leaves_states::Dict{Int, Int}, subtree::Graph{Directed})
  monomial = 1
  
  for edge in edges(subtree)
    dsc = vertex_descendants(subtree, dst(edge))
    dsc_states = filter(pair -> pair.first in dsc, leaves_states)
    group_elem = group_sum(PM, dsc_states)
    monomial = monomial * entry_fourier_parameter(PM, which_group_element(PM, group_elem), edge)
  end
  return monomial
end

function leaves_fourier(PM::GroupBasedPhylogeneticModel{<:PhylogeneticTree}, leaves_states::Dict{Int, Int})
  
  if is_zero_group_sum(PM, leaves_states) 
    return monomial_fourier(PM, leaves_states)
  end 

  S = parameter_ring(PM)[1]
  return S(0)
  
end 

function leaves_fourier(PM::GroupBasedPhylogeneticModel{<:PhylogeneticNetwork}, leaves_states::Dict{Int, Int}, subtree::Graph{Directed})
  
  if is_zero_group_sum(PM, leaves_states) 
    return monomial_fourier(PM, leaves_states, subtree)
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
