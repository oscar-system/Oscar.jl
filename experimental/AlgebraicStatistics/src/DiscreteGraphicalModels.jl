@attributes mutable struct DiscreteGraphicalModel{T, L} <: GraphicalModel{T, L}
  graph::T
  labelings::L
  states::Vector{Int}
  varnames::Dict{Symbol,VarName}
  function DiscreteGraphicalModel(G::Graph{S}, states::Vector{Int}, varnames::Dict{Symbol, VarName}) where S <: GraphTypes
    graph_maps = NamedTuple(_graph_maps(G))
    graph_maps = isempty(graph_maps) ? nothing : graph_maps
    return new{Graph{S}, typeof(graph_maps)}(G, graph_maps, states, varnames)
  end
end

###################################################################################
#
#       Undirected Discrete Graphical Models
#
###################################################################################

@doc raw"""
    discrete_graphical_model(G::Graph{Undirected}, states::Vector{Int}; t_var_name::String="t")

The parametric statistical model associated to an undirected graph.
It contains an undirected graph `G`, a MarkovRing `S` where the vanishing ideal of the model naturally lives, 
and a parameter ring whose variables `t[C](i_1, i_2, ..., i_#C)` correspond to potential functions of each clique. 

# TODO refactor


``` jldoctest 
julia> M = discrete_graphical_model(graph_from_edges([[1,2], [2,3]]), [2,2,2])
Discrete graphical model on an undirected graph with edges
(2, 1), (3, 2) and states [2, 2, 2]

```
"""
function discrete_graphical_model(G::Graph{Undirected}, states::Vector{Int}; t_var_name::String="t", p_var_name::String="p")
  @req length(states) == n_vertices(G) "need one and only one random variable in MarkovRing for each graph vertex"
  DiscreteGraphicalModel(G, states, Dict{Symbol, VarName}(:t => t_var_name, :p => p_var_name))
end

varnames(M::DiscreteGraphicalModel) = M.varnames

states(M::DiscreteGraphicalModel) = M.states

maximal_cliques(M::DiscreteGraphicalModel) = maximal_cliques(M.graph)

function Base.show(io::IO, M::DiscreteGraphicalModel{T, L}) where L
  io = pretty(io)
  if is_terse(io)
    print(io, "Discrete Graphical Model on a")
    !(L == Nothing) && print(io, " labelled")
    print(io, " $T")
  else
    print(io, "Discrete Graphical Model on a $(graph(M))")
  end
  print(io, " with states ", string(M.states))
end

# TODO: Allow the ground field to be specified!
# TODO: Should this return MarkovRing instead which manages the MPolyRing
# together with its generators and any special form they may have.
@doc raw"""
  @attr Tuple{MPolyRing, Vector{QQMPolyRingElem}} model_ring(M::DiscreteGraphicalModel{Graph{Undirected}, T}; cached=false)

Return the ring in which the statistical model lives, together with a
Dict for indexing its generators.

``` jldoctest 
julia> M = discrete_graphical_model(graph_from_edges([[1,2], [2,3]]), [2,2,2])
Discrete graphical model on an undirected graph with edges
(2, 1), (3, 2) and states [2, 2, 2]

julia> model_ring(M)[2]
Dict{Tuple{Int64, Int64, Int64}, QQMPolyRingElem} with 8 entries:
  (1, 2, 1) => p[1, 2, 1]
  (1, 1, 1) => p[1, 1, 1]
  (2, 2, 1) => p[2, 2, 1]
  (1, 2, 2) => p[1, 2, 2]
  (2, 1, 1) => p[2, 1, 1]
  (1, 1, 2) => p[1, 1, 2]
  (2, 2, 2) => p[2, 2, 2]
  (2, 1, 2) => p[2, 1, 2]
```
"""
@attr Tuple{
  MPolyRing,
  GenDict
} function model_ring(M::DiscreteGraphicalModel{Graph{Undirected}, T}; cached=false) where T
  R = markov_ring(M.states...; cached=cached)
  return (ring(R), gens(R))
end

function state_space(M::DiscreteGraphicalModel, C::Vector{Int})
  return length(C) == 1 ? M.states[C[1]] : Iterators.product([1:M.states[i] for i in C]...)
end

@doc raw"""
  @attr Tuple{MPolyRing, GenDict} parameter_ring(M::DiscreteGraphicalModel{Graph{Undirected}, T}; cached=false)

Return the ring of parameters of the statistical model, together with a
Dict for indexing its generators.

``` jldoctest 
julia> M = discrete_graphical_model(graph_from_edges([[1,2], [2,3]]), [2,2,2])
Discrete graphical model on an undirected graph with edges
(2, 1), (3, 2) and states [2, 2, 2]

julia> parameter_ring(M)[2]
Dict{Tuple{Vector{Int64}, Tuple{Int64, Int64}}, QQMPolyRingElem} with 8 entries:
  ([2, 3], (1, 2)) => t[2, 3](1, 2)
  ([2, 1], (1, 2)) => t[2, 1](1, 2)
  ([2, 3], (1, 1)) => t[2, 3](1, 1)
  ([2, 1], (1, 1)) => t[2, 1](1, 1)
  ([2, 3], (2, 2)) => t[2, 3](2, 2)
  ([2, 1], (2, 2)) => t[2, 1](2, 2)
  ([2, 3], (2, 1)) => t[2, 3](2, 1)
  ([2, 1], (2, 1)) => t[2, 1](2, 1)
```
"""
@attr Tuple{
  MPolyRing,
  GenDict
} function parameter_ring(M::DiscreteGraphicalModel{Graph{Undirected}, T}; cached=false) where T
  cliques = sort(sort.(collect.(maximal_cliques(M.graph))))
  Xs = [sort(vec(collect(state_space(M, C)))) for C in cliques]
  params = [(C, x) for (C, X) in Iterators.zip(cliques, Xs) for x in X]
  varnames = sort([M.varnames[:t] * string(C) * string(x) for (C, x) in params])
  R, t = polynomial_ring(QQ, varnames; cached=cached)
  gens_dict = Dict((C, x) => t[i] for (i, (C, x)) in enumerate(params))
  return (R, gens_dict)
end

@doc raw"""
  parameterization(M::DiscreteGraphicalModel{Graph{Undirected}, L})

Creates the polynomial map which parameterizes the vanishing ideal of the
undirected discrete graphical model `M`. It sends each probability generator
`p[s]` of the model ring to the product over all generators `t[C](sp)` of
the parameter ring where `s[C] == sp`, where `C` is a maximal clique of the
graph and `sp` is a marginal state of the random vector.

## Examples

``` jldoctest
julia> M = discrete_graphical_model(graph_from_edges([[1,2], [2,3]]), [2,2,2])
Discrete graphical model on an undirected graph with edges
(2, 1), (3, 2) and states [2, 2, 2]

julia> parametrization(M)
Ring homomorphism
  from multivariate polynomial ring in 8 variables over QQ
  to multivariate polynomial ring in 8 variables over QQ
defined by
  p[1, 1, 1] -> t[1, 2](1, 1)*t[2, 3](1, 1)
  p[2, 1, 1] -> t[1, 2](2, 1)*t[2, 3](1, 1)
  p[1, 2, 1] -> t[1, 2](1, 2)*t[2, 3](2, 1)
  p[2, 2, 1] -> t[1, 2](2, 2)*t[2, 3](2, 1)
  p[1, 1, 2] -> t[1, 2](1, 1)*t[2, 3](1, 2)
  p[2, 1, 2] -> t[1, 2](2, 1)*t[2, 3](1, 2)
  p[1, 2, 2] -> t[1, 2](1, 2)*t[2, 3](2, 2)
  p[2, 2, 2] -> t[1, 2](2, 2)*t[2, 3](2, 2)
```
"""

function parametrization(M::DiscreteGraphicalModel{Graph{Undirected}, L}) where L
  G = graph(M)
  S, pd = model_ring(M)
  R, td = parameter_ring(M)
  # TODO: At this point it would be nice to have the MarkovRing available
  # (e.g., returned from model_ring instead of the inner MPolyRing) so that
  # we could ask what the index s of a given generator "p[s]" is. Right now,
  # we fall back to computing this index.
  images = []
  for p in gens(S)
    s = findfirst(q -> q == p, pd)
    push!(images, prod(td[k] for k in keys(td) if k[2] == s[k[1]]))
  end
  hom(S, R, images)
end

###################################################################################
#
#       Directed Discrete Graphical Models
#
###################################################################################


@doc raw"""
   discrete_graphical_model(G::Graph{Directed}, q_var_name::String="q")

A parametric statistical model associated to a directed acyclic graph.
It contains a directed acylic graph `G`, a MarkovRing `S` where the vanishing ideal of the model naturally lives, 
and a parameter ring whose variables `q[i][j](j_1, j_2, ..., j_#pa(i))` correspond to the conditional distribution 
of the node i given its parents. 

## Examples


``` jldoctest 
julia> M = graphical_model(graph_from_edges(Directed, [[1,3], [2,3], [3, 4]]), markov_ring("1" => 1:2, "2" => 1:2, "3" => 1:2, "4" => 1:2))
discrete graphical model on a directed graph with edges:
(1, 3), (2, 3), (3, 4)
```
"""
function discrete_graphical_model(G::Graph{Directed}, q_var_name::String="q")
  rvs = random_variables(S)
  R, q = polynomial_ring(QQ, vcat([q_var_name * string([i, j]) * string(par_state)
                                   for i in vertices(G)
                                     for j in state_space(S, rvs[i])
                                       for par_state in parental_state_space(G, i, S)], ["h"]))

  param_dict = Dict([[(i, j), Dict()] for i in vertices(G) for j in state_space(S, rvs[i])])
  t = 1
  for i in vertices(G)
    for j in state_space(S, rvs[i])
      for par_state in parental_state_space(G, i, S)
        
        param_dict[i, j][Tuple(par_state)] = q[t]
        t += 1
      end
    end
  end

  return GraphicalModel(G, S, R, param_dict)
end

# returns the joint states of the parents vertex i in a directed graph G
function parental_state_space(G::Graph{Directed}, i::Integer, S::MarkovRing)
  par = inneighbors(G, i)
  par_rvs = random_variables(S)[par]
  
  return collect(state_space(S, par_rvs))
end

@doc raw"""
  parameterization(M::GraphicalModel{Graph{Directed}, MarkovRing})

Creates the polynomial map which parameterizes the vanishing ideal of the directed discrete graphical model `M`.   
The vanishing ideal of the statistical model is the kernel of this map. This ring map is the pull back of the parameterization $\phi_G$
where every joint probability is given by the product of conditional densities of each node given its parents. 
## Examples

``` jldoctest
julia> M = graphical_model(graph_from_edges(Directed, [[1,3], [2,3]]), markov_ring("1" => 1:2, "2" => 1:2, "3" => 1:2))
discrete graphical model on an directed graph with edges:
(1, 3), (2, 3)

julia> parameterization(M)
Ring homomorphism
from multivariate polynomial ring in 8 variables over QQ
to multivariate polynomial ring in 8 variables over QQ
defined by
p[1, 1, 1] -> t[1, 2](1, 1)*t[2, 3](1, 1)
p[2, 1, 1] -> t[1, 2](2, 1)*t[2, 3](1, 1)
p[1, 2, 1] -> t[1, 2](1, 2)*t[2, 3](2, 1)
p[2, 2, 1] -> t[1, 2](2, 2)*t[2, 3](2, 1)
p[1, 1, 2] -> t[1, 2](1, 1)*t[2, 3](1, 2)
p[2, 1, 2] -> t[1, 2](2, 1)*t[2, 3](1, 2)
p[1, 2, 2] -> t[1, 2](1, 2)*t[2, 3](2, 2)
p[2, 2, 2] -> t[1, 2](2, 2)*t[2, 3](2, 2)
"""
# function parametrization(M::GraphicalModel{Graph{Directed}, MarkovRing})
# 
#   G = graph(M)
#   V = vertices(G)
#   S = ring(M)
#   rvs = random_variables(S)
#   R = param_ring(M)
#   q = param_gens(M)
# 
#   images = []
# 
#   for ind in state_space(S)
# 
#       push!(images, prod(map(i-> q[i, ind[i]][Tuple(ind[inneighbors(G, i)])], V)))
#   end
# 
#   lin_constraints  = []
# 
#   for i in V
#     for par_state in parental_state_space(G, i, S)
# 
#       push!(lin_constraints, last(gens(R)) - sum(map(j -> q[i, j][Tuple(par_state)] , state_space(S, rvs[i]))))
#     end
#   end
#   I = ideal(lin_constraints)
# 
#   return hom(ring(S), R, map(f -> reduce(f, gens(I)), images))
# end
