###################################################################################
#
#       Auxiliary graph functions
#
###################################################################################

interior_nodes(graph::Graph) = findall(>=(1), outdegree(graph))

leaves(graph::Graph) = findall(iszero, outdegree(graph))
n_leaves(graph::Graph) = length(leaves(graph))
n_leaves(pt::PhylogeneticTree) = n_leaves(adjacency_tree(pt))
root(pt::PhylogeneticTree) = _root(adjacency_tree(pt))


function roots(graph::Graph)
  n_parents = [length(inneighbors(graph, v)) for v in 1:n_vertices(graph)]
  return findall(iszero, n_parents)
end

root(graph::Graph) = only(roots(graph))


function sort_edges(graph::Graph, sorted_edges::Union{Vector{Edge}, Nothing} = nothing)
  if !isnothing(sorted_edges) 
    if Set(collect(edges(graph))) == Set(sorted_edges)
      return sorted_edges
    else
      @warn("Expected edges $(collect(edges(graph))) and got $sorted_edges")
    end
  end
  edgs = collect(edges(graph))
  leaves_idx = findall(edge -> dst(edge) in leaves(graph), edgs)
  return edgs[vcat(leaves_idx, setdiff(1:length(edgs), leaves_idx))]
end

sort_edges(N::PhylogeneticNetwork{M,L}, sorted_edges::Union{Vector{Edge}, Nothing} = nothing) where {M, L} = sort_edges(graph(N), sorted_edges)
sort_edges(pt::PhylogeneticTree, sorted_edges::Union{Vector{Edge}, Nothing} = nothing) = sort_edges(adjacency_tree(pt), sorted_edges)

descendants(pt::PhylogeneticTree, v::Int) = descendants(adjacency_tree(pt), v)

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

  length(roots) != 1 && return false
  length(lvs) == 0 && return false
  length(vertices(G)) != length(vcat(roots, lvs, tree_nodes, h_nodes)) &&  return false

  return true
end

function hybrid_vertices(G::Graph{Directed})
  A = matrix(ZZ, adjacency_matrix(G))
  indegrees = sum(A[i,:] for i in 1:size(A)[1])
  outdegrees = sum(A[:,i] for i in 1:size(A)[1])

  return findall((indegrees .== 2) .& (outdegrees .== 1))
end

function hybrid_edges(G::Graph{Directed})
  hybrid_nodes = hybrid_vertices(G)
  return [[Edge(j, i) for j in sort(inneighbors(G, i))] for i in hybrid_nodes]
end

function hybrid_edges(G::Graph{Directed}, i)
  hybrid_nodes = hybrid_vertices(G)
  if !(i in hybrid_nodes)
    error("$i is not a hybrid node.")
  end
  return [Edge(j, i) for j in sort(inneighbors(G, i))]
end

function level_phylogenetic_network(G::Graph{Directed})
  if !is_phylogenetic_network(G)
    error("$G is not a phylogenetic network.")
  end

  h_nodes = hybrid_vertices(G)
  bicon_comps = biconnected_components(graph_from_edges(Undirected,edges(G)))

  return maximum([length(intersect(component, h_nodes)) for component in bicon_comps])
end

function tree_edges(N::PhylogeneticNetwork)
  hyb = hybrids(N)
  egdes_N = collect(edges(N))
  return egdes_N[findall(e -> !(dst(e) in collect(keys(hyb))), egdes_N)]
end

###################################################################################
#
#       Auxiliary functions to access parameters
#
###################################################################################

"""
    entry_transition_matrix(PM::PhylogeneticModel, i::Int, j::Int, e::Edge)
    entry_transition_matrix(PM::PhylogeneticModel, i::Int, j::Int, u::Int, v::Int)


Return the `(i, j)`-th entry of the transition matrix associated with edge `e` (or equivalently `Edge(u,v)`) in the phylogenetic model `PM`.

"""
function entry_transition_matrix(PM::PhylogeneticModel{<:AbstractGraph{Directed}, L, <: VarName, T}, i::Int, j::Int, e::Edge) where {L, T}
  tr_mat = transition_matrix(PM)
  return parameter_ring(PM)[2][tr_mat[i,j], e]
end

function entry_transition_matrix(PM::PhylogeneticModel{<:AbstractGraph{Directed}, L, <: VarName, T}, i::Int, j::Int, u::Int, v::Int) where {L, T}
  return entry_transition_matrix(PM, i, j, Edge(u,v))
end

function entry_transition_matrix(PM::PhylogeneticModel{<:AbstractGraph{Directed}, L, <: MPolyRingElem, <: T}, i::Int, j::Int, e::Edge) where {L, T}
  tr_mat = transition_matrix(PM)
  return parameter_ring(PM)[2][e](tr_mat[i,j])
end

function entry_transition_matrix(PM::PhylogeneticModel{<:AbstractGraph{Directed}, L, <: MPolyRingElem, T}, i::Int, j::Int, u::Int, v::Int) where {L, T}
  return entry_transition_matrix(PM, i, j, Edge(u,v))
end

"""
    entry_transition_matrix(PM::GroupBasedPhylogeneticModel, i::Int, j::Int, e::Edge)
    entry_transition_matrix(PM::GroupBasedPhylogeneticModel, i::Int, j::Int, u::Int, v::Int)

Return the `(i, j)`-th entry of the transition matrix associated with edge `e` (or equivalently `Edge(u,v)`) of the underlying `PhylogeneticModel` contained within a `GroupBasedPhylogeneticModel`.
"""
entry_transition_matrix(PM::GroupBasedPhylogeneticModel, i::Int, j::Int, e::Edge) = entry_transition_matrix(phylogenetic_model(PM), i, j, e)
entry_transition_matrix(PM::GroupBasedPhylogeneticModel, i::Int, j::Int, u::Int, v::Int) = entry_transition_matrix(phylogenetic_model(PM), i, j, u, v)

"""
    entry_root_distribution(PM::PhylogeneticModel, i::Int)
    entry_root_distribution(PM::GroupBasedPhylogeneticModel, i::Int)

Return the `i`-th entry of the root distribution parameter for the phylogenetic model associated to `PM`.
"""
entry_root_distribution(PM::PhylogeneticModel, i::Int) = parameter_ring(PM)[3][i]

entry_root_distribution(PM::GroupBasedPhylogeneticModel, i::Int) = entry_root_distribution(phylogenetic_model(PM), i)

"""
    entry_fourier_parameter(PM::GroupBasedPhylogeneticModel, i::Int, e::Edge)

Return the `i`-th Fourier parameter associated with edge `e` (or equivalently `Edge(u,v)`) in a group-based phylogenetic model.
"""
function entry_fourier_parameter(PM::GroupBasedPhylogeneticModel, i::Int, e::Edge)
  x = fourier_parameters(PM)
  return parameter_ring(PM)[2][x[i], e]
end

entry_fourier_parameter(PM::GroupBasedPhylogeneticModel, i::Int, u::Int, v::Int) = entry_fourier_parameter(PM, i, Edge(u,v))

"""
    entry_hybrid_parameter(PM::Union{GroupBasedPhylogeneticModel{<: PhylogeneticNetwork}, PhylogeneticModel{<: PhylogeneticNetwork}}, e::Edge)
    entry_hybrid_parameter(PM::Union{GroupBasedPhylogeneticModel{<: PhylogeneticNetwork}, PhylogeneticModel{<: PhylogeneticNetwork}}, u::Int, v::Int)

Return the parameter associated with a *hybrid edge* `e` (or edge `Edge(u,v)`) in a phylogenetic model or groupd-based model `PM`.

"""
function entry_hybrid_parameter(PM::PhylogeneticModel{<: PhylogeneticNetwork}, e::Edge)
  @req dst(e) in keys(hybrids(graph(PM))) "Edge $e id not a hybrid edge in the Phylogenetic network $(graph(PM))"

  return parameter_ring(PM)[4][e]
end

function entry_hybrid_parameter(PM::GroupBasedPhylogeneticModel{<: PhylogeneticNetwork}, e::Edge) 
  @req dst(e) in keys(hybrids(graph(PM))) "Edge $e id not a hybrid edge in the Phylogenetic network $(graph(PM))"
  return parameter_ring(PM)[3][e]
end

entry_hybrid_parameter(PM::Union{GroupBasedPhylogeneticModel{<: PhylogeneticNetwork}, PhylogeneticModel{<: PhylogeneticNetwork}}, u::Int, v::Int) = entry_hybrid_parameter(PM, Edge(u,v))

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

function fully_observed_probability(PM::PhylogeneticModel, vertices_states::Dict{Int, Int}, tree::AbstractGraph{Directed}) 
  r = root(tree)
  coeff = entry_root_distribution(PM, vertices_states[r])
  exponent = zeros(Int, ngens(parameter_ring(PM)[1]))
  for edge in edges(tree)
    state_parent = vertices_states[src(edge)]
    state_child = vertices_states[dst(edge)]
    gen = entry_transition_matrix(PM, state_parent, state_child, edge)
    exponent += first(exponents(gen))
  end
  return coeff, exponent
end

function leaves_probability(PM::PhylogeneticModel, leaves_states::Dict{Int, Int}, tree::AbstractGraph{Directed})
  int_nodes = interior_nodes(tree)

  interior_indices = collect.(Iterators.product([collect(1:n_states(PM)) for _ in int_nodes]...))  
  vertices_states = leaves_states

  poly = MPolyBuildCtx(parameter_ring(PM)[1])
  for labels in interior_indices
    for (int_node, label) in zip(int_nodes, labels)
      vertices_states[int_node] = label
    end
    push_term!(poly, fully_observed_probability(PM, vertices_states, tree)...)
  end 
  return finish(poly)
end 

function monomial_fourier(PM::GroupBasedPhylogeneticModel{<:PhylogeneticTree}, leaves_states::Dict{Int, Int})
  gr = graph(PM)
  monomial = 1
  
  for edge in edges(gr)
    dsc = descendants(gr, dst(edge))
    dsc_states = filter(pair -> pair.first in dsc, leaves_states)
    group_elem = group_sum(PM, dsc_states)
    monomial = monomial * entry_fourier_parameter(PM, which_group_element(PM, group_elem), edge)
  end
  return monomial
end

function monomial_fourier(PM::GroupBasedPhylogeneticModel{<:PhylogeneticNetwork}, leaves_states::Dict{Int, Int}, subtree::Graph{Directed})
  monomial = 1
  
  for edge in edges(subtree)
    dsc = descendants(subtree, dst(edge))
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
  return zero(S)
end 

function leaves_fourier(PM::GroupBasedPhylogeneticModel{<:PhylogeneticNetwork}, leaves_states::Dict{Int, Int}, subtree::Graph{Directed})
  if is_zero_group_sum(PM, leaves_states) 
    return monomial_fourier(PM, leaves_states, subtree)
  end 

  S = parameter_ring(PM)[1]
  return zero(S)
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

is_zero_group_sum(PM::GroupBasedPhylogeneticModel, states::Dict{Int, Int}) = iszero(group_sum(PM, states))

function which_group_element(PM::GroupBasedPhylogeneticModel, elem::FinGenAbGroupElem)
  G = group(PM)
  return findall([all(g==elem) for g in G])[1]
end
