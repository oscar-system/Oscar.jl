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

function is_group_based_model(matrices::Dict{Edge, MatElem{T}}, G::FinGenAbGroup) where T <: MPolyRingElem
  gb_model = true
  edgs = collect(keys(matrices))
  for e in edgs
    if !is_group_based_matrix(matrices[e], G); gb_model = false; end
  end
  gb_model
end

function is_group_based_matrix(M, G)
  error_occurred = true
  try
    f_group_based_matrix(M, G)
  catch e
    error_occurred = false
  end
  return error_occurred
end

function f_group_based_matrix(M, G)
  ns = ncols(M)
  g_elems = collect(G)
  f = Dict{FinGenAbGroupElem, typeof(matrices[edgs[1]][1,1])}()
  for i in 1:ns
    for j in 1:ns
      if haskey(f, g_elems[i] - g_elems[j]) && f[g_elems[i] - g_elems[j]] != M[i,j]; error(M, " does not have a group structure"); end
      f[g_elems[i] - g_elems[j]] = M[i,j]
    end
  end
  return f
end

function discrete_fourier_transform(M::MatElem, G::FinGenAbGroup; fourier_param_name::String = "x", fourier_param_idx::Int = 0)
  X = character_table(G)
  g_elems = collect(G)
  f = f_group_based_matrix(M, G)

  f̂ = Dict{FinGenAbGroupElem, MPolyRingElem}()
  for i in 1:length(g_elems)
    f̂[g_elems[i]] = sum([QQ(X[i,j])*f[g_elems[j]] for j in 1:length(g_elems)])
  end

  unique_vals = unique(collect(values(f̂)))
  if fourier_param_idx != 0
    S, x = polynomial_ring(QQ, fourier_param_name => (fourier_param_idx:fourier_param_idx, 1:length(unique_vals)); cached=false)
    fourier_param = [x[1, findfirst(unique_vals .== f̂[g])] for g in g_elems]
  else
    S, x = polynomial_ring(QQ, fourier_param_name => (1:length(unique_vals)); cached=false)
    fourier_param = [x[findfirst(unique_vals .== f̂[g])] for g in g_elems]
  end
  return fourier_param
end

function fourier_params_from_matrices(matrices::Dict{Edge, MatElem{T}}, G::FinGenAbGroup; F::Field=QQ) where T <: MPolyRingElem
  edgs = sort_edges(graph_from_edges(Directed,collect(keys(matrices))))
  
  f_params = []
  for i in 1:length(edgs)
    print(i)
    M = matrices[edgs[i]]
    xe = discrete_fourier_transform(M, G; fourier_param_name, fourier_param_idx=i)
    append!(f_params,[xe])
  end

  f_params_list = unique(vcat(f_params...))
  S, x = polynomial_ring(F, ["$var" for var in f_params_list])

  fourier_param = Dict{Edge, Vector{MPolyRingElem}}()
  for i in 1:length(edgs)
    R = base_ring(ideal(f_params[i]))
    inject_xe = hom(R, S, x[[findfirst(y .== S.S) for y in R.S]])
    fourier_param[edgs[i]] = [inject_xe(y) for y in f_params[i]]    
  end

  return(S, fourier_param)
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
  return findall(x -> x == 0, n_parents)[1]
end

########################
#### SORT FUNCTIONS ####
########################
function sort_edges(graph::Graph)
  edgs = collect(edges(graph))
  leaves_idx = findall(edge -> dst(edge) in Oscar.leaves(graph), edgs)
  return edgs[vcat(leaves_idx, setdiff(1:length(edgs), leaves_idx))]
end
