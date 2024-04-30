struct PhylogeneticModel
  graph::Graph{Directed}
  n_states::Int
  prob_ring::MPolyRing{QQFieldElem}
  root_distr::Vector{Any}
  trans_matrices::Dict{Edge, MatElem{QQMPolyRingElem}}
end

function Base.show(io::IO, pm::PhylogeneticModel)
  gr = graph(pm)
  nl = length(leaves(gr))
  ne = length(collect(edges(gr)))
  x = gens(probability_ring(pm))
  print(io, "Phylogenetic model on a tree with $(nl) leaves and $(nl) edges")
  #print(io, " with transition matrices of the form ")
end

struct GroupBasedPhylogeneticModel
  phylo_model::PhylogeneticModel
  fourier_ring::MPolyRing{QQFieldElem}
  fourier_params::Dict{Edge, Vector{QQMPolyRingElem}}
  group::Vector{Vector{Int64}}
end

function Base.show(io::IO, pm::GroupBasedPhylogeneticModel)
  gr = graph(pm)
  nl = length(leaves(gr))
  edgs = collect(edges(gr))
  ne = length(edgs)
  x = gens(probability_ring(pm))
  # M = pm.phylo_model.trans_matrices[edges[1]]
  print(io, "Group-based phylogenetic model on a tree with $(nl) leaves and $(ne) edges")
  #print(io, " with transition matrices of the form ")
end

phylogenetic_model(pm::GroupBasedPhylogeneticModel) =  pm.phylo_model

# calling elements of the structure
@doc raw"""
  graph(pm::PhylogeneticModel)

Return the graph of a PhylogeneticModel `pm`.
"""
graph(pm::PhylogeneticModel) = pm.graph
graph(pm::GroupBasedPhylogeneticModel) = pm.phylo_model.graph

@doc raw"""
  number_states(pm::PhylogeneticModel)

Return the number of states of the PhylogeneticModel `pm`.
"""
number_states(pm::PhylogeneticModel) = pm.n_states
number_states(pm::GroupBasedPhylogeneticModel) = pm.phylo_model.n_states

@doc raw"""
  transition_matrices(pm::PhylogeneticModel)

Returns a dictionary between the edges of the tree specifying a PhylogeneticModel `pm` and attached transition matrices.
"""
transition_matrices(pm::PhylogeneticModel) = pm.trans_matrices
transition_matrices(pm::GroupBasedPhylogeneticModel) = pm.phylo_model.trans_matrices

@doc raw"""
  probability_ring(pm::PhylogeneticModel)

Returns the ring of probability coordinates of the PhylogeneticModel `pm`.
"""
probability_ring(pm::PhylogeneticModel) = pm.prob_ring
probability_ring(pm::GroupBasedPhylogeneticModel) = pm.phylo_model.prob_ring


root_distribution(pm::PhylogeneticModel) = pm.root_distr
root_distribution(pm::GroupBasedPhylogeneticModel) = pm.phylo_model.root_distr


@doc raw"""
  fourier_param(pm::GroupBasedPhylogeneticModel)

Return the Fourier parameters of the GroupBasedPhylogeneticModel `pm` as a vector of eigenvalues of the transition matrices.
"""
fourier_parameters(pm::GroupBasedPhylogeneticModel) = pm.fourier_params

@doc raw"""
  fourier_ring(pm::GroupBasedPhylogeneticModel)

  Returns the ring of Fourier coordinates of the PhylogeneticModel `pm`.
"""
fourier_ring(pm::GroupBasedPhylogeneticModel) = pm.fourier_ring

group_of_model(pm::GroupBasedPhylogeneticModel) = pm.group


#define group-based models

function cavender_farris_neyman_model(graph::Graph{Directed})
  ns = 2
  ne = n_edges(graph)
  R, list_a, list_b = polynomial_ring(QQ, :a => 1:ne, :b => 1:ne; cached=false)
  
  root_distr = repeat([1//ns], outer = ns)
  matrices = Dict{Edge, MatElem}(e => matrix(R, [
    a b 
    b a]) for (a,b,e) in zip(list_a, list_b, edges(graph))
  )
  
  S, list_x = polynomial_ring(QQ, :x => (1:ne, 1:2); cached=false)
  fourier_param = Dict{Edge, Vector{QQMPolyRingElem}}(e => 
    [list_x[i,1], list_x[i,2]] for (i, e) in zip(1:ne, edges(graph)))
  
  group = [[0],[1]]

  pm = PhylogeneticModel(graph, ns, R, root_distr, matrices)
  return GroupBasedPhylogeneticModel(pm, S, fourier_param, group)
end

@doc raw"""
  jukes_cantor_model(graph::Graph{Directed})

Creates a `PhylogeneticModel` based on `graph` whose transition matrices are Jukes Cantor matrices. 

# Examples
```jldoctest
julia> pm = jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]));
```
We can recover the information of `pm`
```jldoctest
julia> graph(pm)
Directed graph with 4 nodes and the following edges:
(4, 1)(4, 2)(4, 3)

julia> number_states(pm)
4

julia> prob_ring(pm)
fourier_ring(pm)
trans_matrices(pm)
fourier_params(pm)
group(pm)
```
"""
function jukes_cantor_model(graph::Graph{Directed})
  ns = 4
  ne = n_edges(graph)
  R, list_a, list_b = polynomial_ring(QQ, :a => 1:ne, :b => 1:ne; cached=false)
  
  root_distr = repeat([1//ns], outer = ns)
  matrices = Dict{Edge, MatElem}(e => matrix(R, [
    a b b b
    b a b b
    b b a b
    b b b a]) for (a,b,e) in zip(list_a, list_b, edges(graph))
  )
  
  S, list_x = polynomial_ring(QQ, :x => (1:ne, 1:2); cached=false)
  fourier_param = Dict{Edge, Vector{QQMPolyRingElem}}(e => 
    [list_x[i,1], list_x[i,2], list_x[i,2], list_x[i,2]] for (i, e) in zip(1:ne, edges(graph)))
  
  group = [[0,0], [0,1], [1,0], [1,1]]

  pm = PhylogeneticModel(graph, ns, R, root_distr, matrices)
  return GroupBasedPhylogeneticModel(pm, S, fourier_param, group)
end

# Kimura 2
function kimura2_model(graph::Graph{Directed})
  ns = 4
  ne = n_edges(graph)
  R, list_a, list_b, list_c = polynomial_ring(QQ, :a => 1:ne, :b => 1:ne, :c => 1:ne; cached=false)
  
  root_distr = repeat([1//ns], outer = ns)
  matrices = Dict{Edge, MatElem}(e => matrix(R, [
    a b c b
    b a b c
    c b a b
    b c b a]) for (a,b,c,e) in zip(list_a, list_b, list_c, edges(graph))
  )

  S, list_x = polynomial_ring(QQ, :x => (1:ne, 1:3); cached=false)
  fourier_param = Dict{Edge, Vector{QQMPolyRingElem}}(e => 
    [list_x[i,1], list_x[i,3], list_x[i,2], list_x[i,2]] for (i, e) in zip(1:ne, edges(graph)))
  
  group = [[0,0], [0,1], [1,0], [1,1]]

  pm = PhylogeneticModel(graph, ns, R, root_distr, matrices)
  return GroupBasedPhylogeneticModel(pm, S, fourier_param, group)
end

# Kimura 3
function kimura3_model(graph::Graph{Directed})
  ns = 4
  ne = n_edges(graph)
  R, list_a, list_b , list_c, list_d= polynomial_ring(QQ, :a => 1:ne, :b => 1:ne, :c => 1:ne, :d => 1:ne; cached=false)
  
  root_distr = repeat([1//ns], outer = ns)
  matrices = Dict{Edge, MatElem}(e => matrix(R, [
    a b c d
    b a d c
    c d a b
    d c b a]) for (a,b,c,d,e) in zip(list_a, list_b, list_c, list_d, edges(graph))
  )

  S, list_x = polynomial_ring(QQ, :x => (1:ne, 1:4); cached=false)
  fourier_param = Dict{Edge, Vector{QQMPolyRingElem}}(e => 
    [list_x[i,1], list_x[i,2], list_x[i,3], list_x[i,4]] for (i, e) in zip(1:ne, edges(graph)))
  
  group = [[0,0], [0,1], [1,0], [1,1]]

  pm = PhylogeneticModel(graph, ns, R, root_distr, matrices)
  return GroupBasedPhylogeneticModel(pm, S, fourier_param, group)
end

# general Markov model
function general_markov_model(graph::Graph{Directed}; number_states = 4)
  ns = number_states
  ne = n_edges(graph)
  R, root_distr, list_m = polynomial_ring(QQ, :π => 1:ns, :m => (1:(ne^2),1:(ns^2)); cached=false)

  matrices = Dict{Edge, MatElem}(e => matrix(R, transpose(reshape(list_m[i,:], ns, ns))) for (i,e) in zip(1:ne, edges(graph)))
  
  return PhylogeneticModel(graph, ns, R, root_distr, matrices)
end
