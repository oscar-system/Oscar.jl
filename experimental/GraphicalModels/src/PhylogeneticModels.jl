struct PhylogeneticModel
  graph::Graph{Directed}
  n_states::Int
  prob_ring::MPolyRing{QQFieldElem}
  root_distr::Vector{Any}
  trans_matrices::Dict{Edge, MatElem{QQMPolyRingElem}}
end

function Base.show(io::IO, pm::PhylogeneticModel)
  gr = graph(pm)
  ns = number_states(pm)
  nl = length(leaves(gr))
  ne = length(collect(edges(gr)))
  root_dist = join(Oscar.root_distribution(pm), ", " )
  M = collect(values(transition_matrices(pm)))[1]
  print(io, "Phylogenetic model on a tree with $(nl) leaves and $(ne) edges \n with distribution at the root [$(root_dist)] \n")
  print(io, " and transition matrix associated to edge i of the form \n ")
  #print_matrix(io, M, ns)
  idx = string(split(string(M[1,1]), "[")[2][1])
  print(io, replace(replace(string(M), "["*idx => "[i"), ";" => ";\n "))
  print(io, ". ")
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
  ne = length(collect(edges(gr)))
  root_dist = join(Oscar.root_distribution(pm), ", " )
  M = collect(values(transition_matrices(pm)))[1]
  idx = string(split(string(M[1,1]), "[")[2][1])
  
  print(io, "Group-based phylogenetic model on a tree with $(nl) leaves and $(ne) edges \n with distribution at the root [$(root_dist)]. \n")
  print(io, " The transition matrix associated to edge i is of the form \n ")
  #print_matrix(io, M, ns)
  print(io, replace(replace(string(M), "["*idx => "[i"), ";" => ";\n "))
  print(io, ", \n and the Fourier parameters are ")
  fp = transpose(collect(values(Oscar.fourier_parameters(pm)))[1])
  fp = replace(string(fp), "QQMPolyRingElem" => "")
  print(io, replace(replace(replace(string(fp), "["*idx => "[i"), ";" => ";\n "), "]]" => "]]."))
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
@doc raw"""
  cavender_farris_neyman_model(graph::Graph{Directed})

Creates a `PhylogeneticModel` based on `graph` whose transition matrices are of type Cavender-Farris-Neyman. 

# Examples
```jldoctest
julia> pm = cavender_farris_neyman_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
Group-based phylogenetic model on a tree with 3 leaves and 3 edges
julia> graph(pm)
Directed graph with 4 nodes and the following edges:
(4, 1)(4, 2)(4, 3)
julia> number_states(pm)
2
julia> probability_ring(pm)
Multivariate polynomial ring in 6 variables a[1], a[2], a[3], b[1], ..., b[3] over rational field
julia> root_distribution(pm)
2-element Vector{Any}:
 1//2
 1//2
julia> transition_matrices(pm)
Dict{Edge, MatElem{QQMPolyRingElem}} with 3 entries:
  Edge(4, 1) => [a[1] b[1]; b[1] a[1]]
  Edge(4, 2) => [a[2] b[2]; b[2] a[2]]
  Edge(4, 3) => [a[3] b[3]; b[3] a[3]]
julia> fourier_ring(pm)
Multivariate polynomial ring in 6 variables x[1, 1], x[2, 1], x[3, 1], x[1, 2], ..., x[3, 2] over rational field
julia> fourier_parameters(pm)
Dict{Edge, Vector{QQMPolyRingElem}} with 3 entries:
  Edge(4, 1) => [x[1, 1], x[1, 2]]
  Edge(4, 2) => [x[2, 1], x[2, 2]]
  Edge(4, 3) => [x[3, 1], x[3, 2]]
julia> group_model(pm)
???
```
"""
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
julia> pm = jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
Group-based phylogenetic model on a tree with 3 leaves and 3 edges
julia> graph(pm)
Directed graph with 4 nodes and the following edges:
(4, 1)(4, 2)(4, 3)
julia> number_states(pm)
4
julia> probability_ring(pm)
Multivariate polynomial ring in 6 variables a[1], a[2], a[3], b[1], ..., b[3] over rational field
julia> root_distribution(pm)
4-element Vector{Any}:
 1//4
 1//4
 1//4
 1//4
julia> transition_matrices(pm)
Dict{Edge, MatElem{QQMPolyRingElem}} with 3 entries:
  Edge(4, 1) => [a[1] b[1] b[1] b[1]; b[1] a[1] b[1] b[1]; b[1] b[1] a[1] b[1]; b[1] b[1] b[1] a[1]]
  Edge(4, 2) => [a[2] b[2] b[2] b[2]; b[2] a[2] b[2] b[2]; b[2] b[2] a[2] b[2]; b[2] b[2] b[2] a[2]]
  Edge(4, 3) => [a[3] b[3] b[3] b[3]; b[3] a[3] b[3] b[3]; b[3] b[3] a[3] b[3]; b[3] b[3] b[3] a[3]]
julia> fourier_ring(pm)
Multivariate polynomial ring in 6 variables x[1, 1], x[2, 1], x[3, 1], x[1, 2], ..., x[3, 2] over rational field
julia> fourier_parameters(pm)
Dict{Edge, Vector{QQMPolyRingElem}} with 3 entries:
  Edge(4, 1) => [x[1, 1], x[1, 2], x[1, 2], x[1, 2]]
  Edge(4, 2) => [x[2, 1], x[2, 2], x[2, 2], x[2, 2]]
  Edge(4, 3) => [x[3, 1], x[3, 2], x[3, 2], x[3, 2]]
julia> group_model(pm)
???
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
@doc raw"""
    kimura2_model(graph::Graph{Directed})

Creates a `PhylogeneticModel` based on `graph` whose transition matrices are Kimura 2-parameter matrices. 

# Examples
```jldoctest
julia> pm = kimura2_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
Group-based phylogenetic model on a tree with 3 leaves and 3 edges
julia> graph(pm)
Directed graph with 4 nodes and the following edges:
(4, 1)(4, 2)(4, 3)
julia> number_states(pm)
4
julia> probability_ring(pm)
Multivariate polynomial ring in 9 variables a[1], a[2], a[3], b[1], ..., c[3] over rational field
julia> root_distribution(pm)
4-element Vector{Any}:
 1//4
 1//4
 1//4
 1//4
julia> transition_matrices(pm)
Dict{Edge, MatElem{QQMPolyRingElem}} with 3 entries:
  Edge(4, 1) => [a[1] b[1] c[1] b[1]; b[1] a[1] b[1] c[1]; c[1] b[1] a[1] b[1]; b[1] c[1] b[1] a[1]]
  Edge(4, 2) => [a[2] b[2] c[2] b[2]; b[2] a[2] b[2] c[2]; c[2] b[2] a[2] b[2]; b[2] c[2] b[2] a[2]]
  Edge(4, 3) => [a[3] b[3] c[3] b[3]; b[3] a[3] b[3] c[3]; c[3] b[3] a[3] b[3]; b[3] c[3] b[3] a[3]]
julia> fourier_ring(pm)
Multivariate polynomial ring in 9 variables x[1, 1], x[2, 1], x[3, 1], x[1, 2], ..., x[3, 3] over rational field
julia> fourier_parameters(pm)
Dict{Edge, Vector{QQMPolyRingElem}} with 3 entries:
  Edge(4, 1) => [x[1, 1], x[1, 2], x[1, 2], x[1, 2]]
  Edge(4, 2) => [x[2, 1], x[2, 2], x[2, 2], x[2, 2]]
  Edge(4, 3) => [x[3, 1], x[3, 2], x[3, 2], x[3, 2]]
julia> group_model(pm)
???
```
"""
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
@doc raw"""
    kimura3_model(graph::Graph{Directed})

Creates a `PhylogeneticModel` based on `graph` whose transition matrices are Kimura 3-parameter matrices. 

# Examples
```jldoctest
julia> pm = kimura3_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
Group-based phylogenetic model on a tree with 3 leaves and 3 edges
julia> graph(pm)
Directed graph with 4 nodes and the following edges:
(4, 1)(4, 2)(4, 3)
julia> number_states(pm)
4
julia> probability_ring(pm)
Multivariate polynomial ring in 12 variables a[1], a[2], a[3], b[1], ..., d[3] over rational field
julia> root_distribution(pm)
4-element Vector{Any}:
 1//4
 1//4
 1//4
 1//4
julia> transition_matrices(pm)
Dict{Edge, MatElem{QQMPolyRingElem}} with 3 entries:
  Edge(4, 1) => [a[1] b[1] c[1] d[1]; b[1] a[1] d[1] c[1]; c[1] d[1] a[1] b[1]; d[1] c[1] b[1] a[1]]
  Edge(4, 2) => [a[2] b[2] c[2] d[2]; b[2] a[2] d[2] c[2]; c[2] d[2] a[2] b[2]; d[2] c[2] b[2] a[2]]
  Edge(4, 3) => [a[3] b[3] c[3] d[3]; b[3] a[3] d[3] c[3]; c[3] d[3] a[3] b[3]; d[3] c[3] b[3] a[3]]
julia> fourier_ring(pm)
Multivariate polynomial ring in 12 variables x[1, 1], x[2, 1], x[3, 1], x[1, 2], ..., x[3, 4] over rational field
julia> fourier_parameters(pm)
Dict{Edge, Vector{QQMPolyRingElem}} with 3 entries:
  Edge(4, 1) => [x[1, 1], x[1, 2], x[1, 3], x[1, 4]]
  Edge(4, 2) => [x[2, 1], x[2, 2], x[2, 3], x[2, 4]]
  Edge(4, 3) => [x[3, 1], x[3, 2], x[3, 3], x[3, 4]]
julia> group_model(pm)
???
```
"""
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
@doc raw"""
    general_markov_model(graph::Graph{Directed})

Creates a `PhylogeneticModel` based on `graph` whose transition matrices are stochastic with no further constraints. 

# Examples
```jldoctest
julia> pm = general_markov_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
Phylogenetic model on a tree with 3 leaves and 3 edges
julia> graph(pm)
Directed graph with 4 nodes and the following edges:
(4, 1)(4, 2)(4, 3)
julia> number_states(pm)
4
julia> probability_ring(pm)
Multivariate polynomial ring in 148 variables π[1], π[2], π[3], π[4], ..., m[9, 16] over rational field
julia> root_distribution(pm)
4-element Vector{Any}:
 π[1]
 π[2]
 π[3]
 π[4]
julia> transition_matrices(pm)
Dict{Edge, MatElem{QQMPolyRingElem}} with 3 entries:
  Edge(4, 1) => [m[1, 1] m[1, 2] m[1, 3] m[1, 4]; m[1, 5] m[1, 6] m[1, 7] m[1, 8]; m[1, 9] m[1, 10] m[1, 11] m[1, 12]; m[1, 13] m[1, 14] m[1, 15] m[1, 16]]
  Edge(4, 2) => [m[2, 1] m[2, 2] m[2, 3] m[2, 4]; m[2, 5] m[2, 6] m[2, 7] m[2, 8]; m[2, 9] m[2, 10] m[2, 11] m[2, 12]; m[2, 13] m[2, 14] m[2, 15] m[2, 16]]
  Edge(4, 3) => [m[3, 1] m[3, 2] m[3, 3] m[3, 4]; m[3, 5] m[3, 6] m[3, 7] m[3, 8]; m[3, 9] m[3, 10] m[3, 11] m[3, 12]; m[3, 13] m[3, 14] m[3, 15] m[3, 16]]
```
"""
function general_markov_model(graph::Graph{Directed}; number_states = 4)
  ns = number_states
  ne = n_edges(graph)
  R, root_distr, list_m = polynomial_ring(QQ, :π => 1:ns, :m => (1:(ne^2),1:(ns),1:(ns)); cached=false)

  matrices = Dict{Edge, MatElem}(e => matrix(R, reshape(list_m[i,:,:], ns, ns)) for (i,e) in zip(1:ne, edges(graph)))

  return PhylogeneticModel(graph, ns, R, root_distr, matrices)
end


function affine_phylogenetic_model(pm::PhylogeneticModel)
  gr = graph(pm)
  ns = number_states(pm)
  trans_mat = transition_matrices(pm)
  for e in edges(gr)
    [trans_mat[e][i,i] = 1 - (sum(trans_mat[e][i,:]) - trans_mat[e][i,i]) for i in 1:ns]
  end
  r = root_distribution(pm)
  r[1] = 1 - sum(r[2:ns])
  return PhylogeneticModel(gr, ns, probability_ring(pm), root_distribution(pm), trans_mat)
end

function affine_phylogenetic_model(pm::GroupBasedPhylogeneticModel)
  gr = graph(pm)
  affine_pm = affine_phylogenetic_model(pm.phylo_model)
  fourier_param = fourier_parameters(pm)
  S = fourier_ring(pm)
  for e in edges(gr)
    fourier_param[e][1] = S(1)
  end
  return GroupBasedPhylogeneticModel(affine_pm, S, fourier_param, group_of_model(pm))
end
    
