######################################
#### PHYLOGENETIC DATA STRUCTURES ####
######################################

@attributes mutable struct PhylogeneticModel{L} <: GraphicalModel{Directed, L}
  # do need to for T to be directed here?
  base_field::Field
  graph::Graph{Directed}
  labelings::L
  trans_mat_structure::Matrix{<: VarName}
  model_parameter_name::VarName
  n_states::Int
  root_distribution::Vector
  
  function PhylogeneticModel(F::Field,
                             G::Graph{Directed},
                             trans_mat_structure::Matrix{<: VarName},
                             n_states::Union{Nothing, Int} = nothing,
                             root_distribution::Union{Nothing, Vector} = nothing,
                             varname::VarName="p")
    if isnothing(n_states)
      n_states = 4
    end
    if isnothing(root_distribution)
      root_distribution = repeat([1//n_states], outer = n_states)
    end
    graph_maps = NamedTuple(_graph_maps(G))
    graph_maps = isempty(graph_maps) ? nothing : graph_maps
    return new{typeof(graph_maps)}(F, G,
                                   graph_maps,
                                   trans_mat_structure,
                                   varname,
                                   n_states,
                                   root_distribution)
  end

  function PhylogeneticModel(G::Graph{Directed},
                             trans_mat_structure::Matrix{<: VarName},
                             n_states::Union{Nothing, Int} = nothing,
                             root_distribution::Union{Nothing, Vector} = nothing,
                             varname::VarName="p")
    return PhylogeneticModel(QQ, G, trans_mat_structure, n_states, root_distribution, varname)
  end
end

n_states(PM::PhylogeneticModel) = PM.n_states
trans_matrix(PM::PhylogeneticModel) = PM.trans_mat_structure
base_field(PM::PhylogeneticModel) = PM.base_field
varname(PM::PhylogeneticModel) = PM.model_parameter_name

@attr Tuple{MPolyRing, GenDict} function parameter_ring(PM::PhylogeneticModel; cached=false)
  vars = unique(trans_matrix(PM))
  edge_gens = [x => 1:n_edges(graph(PM)) for x in vars]
  R, x... = polynomial_ring(base_field(PM), edge_gens...; cached=cached)

  R, Dict{Tuple{VarName, Edge}, MPolyRingElem}(
    (vars[i], e) => x[i][j] for i in 1:length(vars), (j,e) in enumerate(edges(graph(PM)))
  )
end

@attr Tuple{
  MPolyRing,
  Array
} function model_ring(PM::PhylogeneticModel; cached=false)
  leave_nodes = leaves(graph(PM))
  leaves_indices = collect(Iterators.product(
    [tuple(1:n_states(PM)...) for _ in leave_nodes]...))

  return polynomial_ring(base_field(PM),
                         ["$(varname(PM))[$(join(x, ", "))]" for x in leaves_indices];
                         cached=cached)
end

function parametrization(PM::PhylogeneticModel)
  g =  graph(PM)
  int_nodes = interior_nodes(g)
  lvs_nodes = leaves(g)
  n_states = n_states(pm)

  _, e_gens = parameter_ring(PM)
  tr_mat = trans_matrix(pm)
  r = root(g)
  monomial = root_dist[states[r]]

  for e in edges(g)
    
  end
end

# prob_ring::MPolyRing{QQFieldElem}
# root_distr::Vector{Any} #this need to become more precise
# trans_matrices::Dict{Edge, MatElem{QQMPolyRingElem}}


function Base.show(io::IO, pm::PhylogeneticModel)
  gr = graph(pm)
  ns = n_states(pm)
  nl = length(leaves(gr))
  ne = length(collect(edges(gr)))
  #root_dist = join(Oscar.root_distribution(pm), ", " )
  #M = collect(values(transition_matrices(pm)))[1]

  # commenting out root distribution for now
  print(io, "Phylogenetic model on a tree with $(nl) leaves and $(ne) edges") # \n with distribution at the root [$(root_dist)] \n")
  print(io, " and transition matrices of the form \n ")

  # printing this matrix can probably be improved
  print(io, "$(pm.trans_mat_structure)")
end

struct GroupBasedPhylogeneticModel
  phylo_model::PhylogeneticModel
  fourier_ring::MPolyRing{QQFieldElem}
  fourier_params::Dict{Edge, Vector{QQMPolyRingElem}}
  group::Vector{FinGenAbGroupElem}
end

function Base.show(io::IO, pm::GroupBasedPhylogeneticModel)
  gr = graph(pm)
  nl = length(leaves(gr))
  ne = length(collect(edges(gr)))
  root_dist = join(Oscar.root_distribution(pm), ", ")
  c_edg = 2
  p_edg = inneighbors(gr, c_edg)[1]
  findall(x-> x==2, dst.(edges(gr)))
  M = transition_matrices(pm)[Edge(p_edg, c_edg)]
  idx = string(split(string(M[1,1]), "[")[2][1])

  print(io, "Group-based phylogenetic model on a tree with $(nl) leaves and $(ne) edges \n with distribution at the root [$(root_dist)]. \n")
  print(io, " The transition matrix associated to edge i is of the form \n ")
  print(io, replace(replace(string(M), "["*idx => "[i"), ";" => ";\n "))
  print(io, ", \n and the Fourier parameters are ")
  fp = transpose(fourier_parameters(pm)[Edge(p_edg, c_edg)])
  fp = replace(string(fp), "QQMPolyRingElem" => "")
  print(io, replace(replace(replace(string(fp), "["*idx => "[i"), ";" => ";\n "), "]]" => "]]."))
end

#########################################################################
#### ATTRIBUTES OF PhylogeneticModel AND GroupBasedPhylogeneticModel ####
#########################################################################
@doc raw"""
    phylogenetic_model(pm::GroupBasedPhylogeneticModel)

Return the complete information of a `PhylogeneticModel` or `GroupBasedPhylogeneticModel` `pm`.

# Examples
```jldoctest
julia> pm = jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]));

julia> phylogenetic_model(pm)
Phylogenetic model on a tree with 3 leaves and 3 edges 
 with distribution at the root [1//4, 1//4, 1//4, 1//4] 
 and transition matrix associated to edge i of the form 
 [a[i] b[i] b[i] b[i];
  b[i] a[i] b[i] b[i];
  b[i] b[i] a[i] b[i];
  b[i] b[i] b[i] a[i]]. 
```
"""
phylogenetic_model(pm::GroupBasedPhylogeneticModel) = pm.phylo_model

@doc raw"""
    graph(pm::PhylogeneticModel)

Return the graph of a `PhylogeneticModel` `pm`.

# Examples
```jldoctest
julia> pm = jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]));

julia> graph(pm)
Directed graph with 4 nodes and the following edges:
(4, 1)(4, 2)(4, 3)
```
"""
graph(pm::PhylogeneticModel) = pm.graph
graph(pm::GroupBasedPhylogeneticModel) = pm.phylo_model.graph

@doc raw"""
    number_states(pm::PhylogeneticModel)

Return the number of states of the `PhylogeneticModel` `pm`.

# Examples
```jldoctest
julia> pm = jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]));

julia> number_states(pm)
4
```
"""
number_states(pm::PhylogeneticModel) = pm.n_states
number_states(pm::GroupBasedPhylogeneticModel) = pm.phylo_model.n_states

@doc raw"""
    transition_matrices(pm::PhylogeneticModel)

Return a dictionary between the edges of the tree specifying the `PhylogeneticModel` `pm` and their attached transition matrices.

# Examples
```jldoctest; filter = r".*"
julia> pm = jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]));

julia> transition_matrices(pm)
Dict{Edge, MatElem{QQMPolyRingElem}} with 3 entries:
  Edge(4, 1) => [a[1] b[1] b[1] b[1]; b[1] a[1] b[1] b[1]; b[1] b[1] a[1] b[1];…
  Edge(4, 2) => [a[2] b[2] b[2] b[2]; b[2] a[2] b[2] b[2]; b[2] b[2] a[2] b[2];…
  Edge(4, 3) => [a[3] b[3] b[3] b[3]; b[3] a[3] b[3] b[3]; b[3] b[3] a[3] b[3];…
```
"""
transition_matrices(pm::PhylogeneticModel) = pm.trans_matrices
transition_matrices(pm::GroupBasedPhylogeneticModel) = pm.phylo_model.trans_matrices

@doc raw"""
    probability_ring(pm::PhylogeneticModel)

Return the ring of probability coordinates of the `PhylogeneticModel` `pm`.

# Examples
```jldoctest
julia> pm = jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]));

julia> probability_ring(pm)
Multivariate polynomial ring in 6 variables a[1], a[2], a[3], b[1], ..., b[3]
  over rational field
```
"""
probability_ring(pm::PhylogeneticModel) = pm.prob_ring
probability_ring(pm::GroupBasedPhylogeneticModel) = pm.phylo_model.prob_ring

@doc raw"""
    root_distribution(pm::PhylogeneticModel)

Return the distribution of the random variable at the root of the tree specifying the `PhylogeneticModel` `pm`.

# Examples
```jldoctest
julia> pm = jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]));

julia> root_distribution(pm)
4-element Vector{Any}:
 1//4
 1//4
 1//4
 1//4
```
"""
root_distribution(pm::PhylogeneticModel) = pm.root_distribution
root_distribution(pm::GroupBasedPhylogeneticModel) = pm.phylo_model.root_distribution

@doc raw"""
    fourier_parameters(pm::GroupBasedPhylogeneticModel)

Return the Fourier parameters of the `GroupBasedPhylogeneticModel` `pm` as a vector of eigenvalues of the transition matrices.

# Examples
```jldoctest; filter = r".*"
julia> pm = jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]));

julia> fourier_parameters(pm)
Dict{Edge, Vector{QQMPolyRingElem}} with 3 entries:
  Edge(4, 1) => [x[1, 1], x[1, 2], x[1, 2], x[1, 2]]
  Edge(4, 2) => [x[2, 1], x[2, 2], x[2, 2], x[2, 2]]
  Edge(4, 3) => [x[3, 1], x[3, 2], x[3, 2], x[3, 2]]
```
"""
fourier_parameters(pm::GroupBasedPhylogeneticModel) = pm.fourier_params

@doc raw"""
    fourier_ring(pm::GroupBasedPhylogeneticModel)

Return the ring of Fourier coordinates of the `PhylogeneticModel` `pm`.

# Examples
```jldoctest
julia> pm = jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]));

julia> fourier_ring(pm)
Multivariate polynomial ring in 6 variables x[1, 1], x[2, 1], x[3, 1], x[1, 2], ..., x[3, 2]
  over rational field
```
"""
fourier_ring(pm::GroupBasedPhylogeneticModel) = pm.fourier_ring

@doc raw"""
    group_of_model(pm::GroupBasedPhylogeneticModel)

Return the group the `GroupBasedPhylogeneticModel` `pm` is based on.

# Examples
```jldoctest
julia> pm = jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]));

julia> group_of_model(pm)
4-element Vector{FinGenAbGroupElem}:
 [0, 0]
 [0, 1]
 [1, 0]
 [1, 1]
```
"""
group_of_model(pm::GroupBasedPhylogeneticModel) = pm.group


############################
#### GROUP-BASED MODELS ####
############################

@doc raw"""
    cavender_farris_neyman_model(graph::Graph{Directed})

Create a `PhylogeneticModel` based on `graph` whose transition matrices are of type Cavender-Farris-Neyman. 

# Examples
```jldoctest
julia> cavender_farris_neyman_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
Group-based phylogenetic model on a tree with 3 leaves and 3 edges 
 with distribution at the root [1//2, 1//2]. 
 The transition matrix associated to edge i is of the form 
 [a[i] b[i];
  b[i] a[i]], 
 and the Fourier parameters are [x[i, 1] x[i, 2]].
```
"""
function cavender_farris_neyman_model(graph::Graph{Directed})
  ns = 2
  ne = n_edges(graph)
  R, list_a, list_b = polynomial_ring(QQ, :a => 1:ne, :b => 1:ne; cached=false)
  
  root_distr = repeat([1//ns], outer = ns)
  edgs = sort_edges(graph)
  matrices = Dict{Edge, MatElem}(e => matrix(R, [
    a b 
    b a]) for (a,b,e) in zip(list_a, list_b, edgs))

  S, list_x = polynomial_ring(QQ, :x => (1:ne, 1:2); cached=false)
  fourier_param = Dict{Edge, Vector{QQMPolyRingElem}}(e => 
    [list_x[i,1], list_x[i,2]] for (i, e) in zip(1:ne, edgs))
  
  G = collect(abelian_group(2))
  group = [G[1],G[2]]

  pm = PhylogeneticModel(graph, ns, R, root_distr, matrices)
  return GroupBasedPhylogeneticModel(pm, S, fourier_param, group)
end

@doc raw"""
    jukes_cantor_model(graph::Graph{Directed})

Create a `PhylogeneticModel` based on `graph` whose transition matrices are Jukes Cantor matrices. 

# Examples
```jldoctest
julia> jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
Group-based phylogenetic model on a tree with 3 leaves and 3 edges 
 with distribution at the root [1//4, 1//4, 1//4, 1//4]. 
 The transition matrix associated to edge i is of the form 
 [a[i] b[i] b[i] b[i];
  b[i] a[i] b[i] b[i];
  b[i] b[i] a[i] b[i];
  b[i] b[i] b[i] a[i]], 
 and the Fourier parameters are [x[i, 1] x[i, 2] x[i, 2] x[i, 2]].
```
"""
function jukes_cantor_model(G::Graph{Directed})
  PhylogeneticModel(G, [:a :b :b :b;
                        :b :a :b :b;
                        :b :b :a :b;
                        :b :b :b :a]
                    )
  #G = collect(abelian_group(2,2))
  #group = [G[1],G[3],G[2],G[4]]
  #
  #pm = PhylogeneticModel(graph, ns, R, root_distr, matrices)
  #return GroupBasedPhylogeneticModel(pm, S, fourier_param, group)
end

@doc raw"""
    kimura2_model(graph::Graph{Directed})

Create a `PhylogeneticModel` based on `graph` whose transition matrices are Kimura 2-parameter matrices. 

# Examples
```jldoctest
julia> kimura2_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
Group-based phylogenetic model on a tree with 3 leaves and 3 edges
 with distribution at the root [1//4, 1//4, 1//4, 1//4]. 
 The transition matrix associated to edge i is of the form 
 [a[i] b[i] c[i] b[i];
  b[i] a[i] b[i] c[i];
  c[i] b[i] a[i] b[i];
  b[i] c[i] b[i] a[i]], 
 and the Fourier parameters are [x[i, 1] x[i, 3] x[i, 2] x[i, 2]].
```
"""
function kimura2_model(graph::Graph{Directed})
  ns = 4
  ne = n_edges(graph)
  R, list_a, list_b, list_c = polynomial_ring(QQ, :a => 1:ne, :b => 1:ne, :c => 1:ne; cached=false)
  
  root_distr = repeat([1//ns], outer = ns)
  edgs = sort_edges(graph)
  matrices = Dict{Edge, MatElem}(e => matrix(R, [
    a b c b
    b a b c
    c b a b
    b c b a]) for (a,b,c,e) in zip(list_a, list_b, list_c, edgs))

  S, list_x = polynomial_ring(QQ, :x => (1:ne, 1:3); cached=false)
  fourier_param = Dict{Edge, Vector{QQMPolyRingElem}}(e => 
    [list_x[i,1], list_x[i,3], list_x[i,2], list_x[i,2]] for (i, e) in zip(1:ne, edgs))
  
  G = collect(abelian_group(2,2))
  group = [G[1],G[3],G[2],G[4]]

  pm = PhylogeneticModel(graph, ns, R, root_distr, matrices)
  return GroupBasedPhylogeneticModel(pm, S, fourier_param, group)
end

@doc raw"""
    kimura3_model(graph::Graph{Directed})

Create a `PhylogeneticModel` based on `graph` whose transition matrices are Kimura 3-parameter matrices. 

# Examples
```jldoctest
julia> kimura3_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
Group-based phylogenetic model on a tree with 3 leaves and 3 edges 
 with distribution at the root [1//4, 1//4, 1//4, 1//4]. 
 The transition matrix associated to edge i is of the form 
 [a[i] b[i] c[i] d[i];
  b[i] a[i] d[i] c[i];
  c[i] d[i] a[i] b[i];
  d[i] c[i] b[i] a[i]], 
 and the Fourier parameters are [x[i, 1] x[i, 2] x[i, 3] x[i, 4]].
```
"""
function kimura3_model(graph::Graph{Directed})
  ns = 4
  ne = n_edges(graph)
  R, list_a, list_b , list_c, list_d= polynomial_ring(QQ, :a => 1:ne, :b => 1:ne, :c => 1:ne, :d => 1:ne; cached=false)
  
  root_distr = repeat([1//ns], outer = ns)
  edgs = sort_edges(graph)
  matrices = Dict{Edge, MatElem}(e => matrix(R, [
    a b c d
    b a d c
    c d a b
    d c b a]) for (a,b,c,d,e) in zip(list_a, list_b, list_c, list_d, edgs))

  S, list_x = polynomial_ring(QQ, :x => (1:ne, 1:4); cached=false)
  fourier_param = Dict{Edge, Vector{QQMPolyRingElem}}(e => 
    [list_x[i,1], list_x[i,2], list_x[i,3], list_x[i,4]] for (i, e) in zip(1:ne, edgs))
  
  G = collect(abelian_group(2,2))
  group = [G[1],G[3],G[2],G[4]]

  pm = PhylogeneticModel(graph, ns, R, root_distr, matrices)
  return GroupBasedPhylogeneticModel(pm, S, fourier_param, group)
end

##############################
#### GENERAL MARKOV MODEL ####
##############################

@doc raw"""
    general_markov_model(graph::Graph{Directed})

Create a `PhylogeneticModel` based on `graph` whose transition matrices are stochastic with no further constraints. 

# Examples
```jldoctest
julia> general_markov_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
Phylogenetic model on a tree with 3 leaves and 3 edges 
 with distribution at the root [π[1], π[2], π[3], π[4]] 
 and transition matrix associated to edge i of the form 
 [m[i, 1, 1] m[i, 1, 2] m[i, 1, 3] m[i, 1, 4];
  m[i, 2, 1] m[i, 2, 2] m[i, 2, 3] m[i, 2, 4];
  m[i, 3, 1] m[i, 3, 2] m[i, 3, 3] m[i, 3, 4];
  m[i, 4, 1] m[i, 4, 2] m[i, 4, 3] m[i, 4, 4]].
```
"""
function general_markov_model(graph::Graph{Directed}; number_states = 4)
  ns = number_states
  ne = n_edges(graph)
  R, root_distr, list_m = polynomial_ring(QQ, :π => 1:ns, :m => (1:(ne^2),1:(ns),1:(ns)); cached=false)

  edgs = sort_edges(graph)
  matrices = Dict{Edge, MatElem}(e => matrix(R, reshape(list_m[i,:,:], ns, ns)) for (i,e) in zip(1:ne, edgs))

  return PhylogeneticModel(graph, ns, R, root_distr, matrices)
end


##############################
#### AFFINE PHYLO MODELS #####
##############################

@doc raw"""
    affine_phylogenetic_model!(pm::PhylogeneticModel)

Moves a `PhylogeneticModel` or `GroupBasedPhylogeneticModel` from projective into affine space.

# Examples
```jldoctest
julia> pm = jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]));

julia> affine_phylogenetic_model!(pm)
Group-based phylogenetic model on a tree with 3 leaves and 3 edges
 with distribution at the root [1//4, 1//4, 1//4, 1//4].
 The transition matrix associated to edge i is of the form
 [-3*b[i]+1 b[i] b[i] b[i];
  b[i] -3*b[i]+1 b[i] b[i];
  b[i] b[i] -3*b[i]+1 b[i];
  b[i] b[i] b[i] -3*b[i]+1],
 and the Fourier parameters are [1 x[i, 2] x[i, 2] x[i, 2]].
```
"""
function affine_phylogenetic_model!(pm::PhylogeneticModel)
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

function affine_phylogenetic_model!(pm::GroupBasedPhylogeneticModel)
  gr = graph(pm)
  affine_pm = affine_phylogenetic_model!(pm.phylo_model)
  fourier_param = fourier_parameters(pm)
  S = fourier_ring(pm)
  for e in edges(gr)
    fourier_param[e][1] = S(1)
  end
  return GroupBasedPhylogeneticModel(affine_pm, S, fourier_param, group_of_model(pm))
end

#######################################################
####  PHYLOGENETIC MODEL FROM TRANSITION MATRICES #####
#######################################################

function phylogenetic_model_from_matrices(matrices::Dict{Edge, MatElem{T}};
                                          root_distr::Vector=[]) where T <: MPolyRingElem
  phylogenetic_model_input_check(matrices; root_distr=root_distr)
  edgs = collect(keys(matrices))
  ns = ncols(matrices[edgs[1]])
  F = coefficient_ring(base_ring(matrices[edgs[1]]))
  vars_mat = vec(unique([string(x) for x in hcat(hcat(collect(values(matrices))...)...)]))
  if length(root_distr) == 0
    R, root_distr, = polynomial_ring(F, :π => 1:ns, vars_mat) # this :π should be changed, to :pi?
  elseif isa(root_distr[1], Number)
    R, = polynomial_ring(F, vars_mat)
  else
    R, = polynomial_ring(F, [string(x) for x in root_distr], vars_mat)
  end

  PhylogeneticModel(graph_from_edges(Directed,edgs), ns, R, root_distr, matrices)
end

function phylogenetic_model_input_check(matrices::Dict{Edge, MatElem{T}};
                                        root_distr::Vector=[]) where T <: MPolyRingElem
  edgs = collect(keys(matrices))
  ns = unique([ncols(matrices[edgs[i]]) for i in 1:length(edgs)])

  @req !all(is_square.(collect(values(matrices)))) "Not all matrices are square"
  @req length(ns) != 1 "Matrices of different size"
  @req length(unique(typeof.(collect(values(matrices))))) > 1 "Matrices of different type"
  @req length(unique(base_ring.(collect(values(matrices))))) > 1 "Matrices defined in different polynomial rings"
  @req !isempty(root_distr) && length(root_distr) != ns[1] "Different number of states on transition matrices and and distribution at the root"
  @req length(unique(typeof.(root_distr))) > 1 "Entries of root distribution of different type"
  if !isempty(root_distr)
    r = root_distr[1]; m = matrices[edgs[1]][1,1]
    if !isa(r, Number)
      F_m = coefficient_ring(m)
      F_r = coefficient_ring(r)
      @req F_m != F_r "Transition matrices and distribution at the root defined in different fields"
    end
  end
end

function group_based_phylogenetic_model_from_matrices(matrices::Dict{Edge, MatElem{T}}, G::FinGenAbGroup;
                                                      root_distr::Vector=[], fourier_param_name::String = "x") where T <: MPolyRingElem
  # fourier_param_name should probably be changed to a VarName?
  # phylogenetic_model_input_check(matrices; root_distr=root_distr)
  is_group_based_model(matrices, G)
  edgs = collect(keys(matrices))
  ns = ncols(matrices[edgs[i]]) 
  if length(root_distr) == 0
    pm = phylogenetic_model_from_matrices(matrices; root_distr=repeat([1 / ns], ns))
  else
    pm = phylogenetic_model_from_matrices(matrices; root_distr=root_distr)
  end
  (S, fourier_param) = fourier_params_from_matrices(matrices, G)
  GroupBasedPhylogeneticModel(pm, S, fourier_param, collect(G))
end

