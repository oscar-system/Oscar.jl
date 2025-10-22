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

``` jldoctest 
julia> M = discrete_graphical_model(graph_from_edges([[1,2], [2,3]]), [2,2,2])
Discrete Graphical Model on a Undirected graph with 3 nodes and 2 edges with states [2, 2, 2]
```
"""
function discrete_graphical_model(G::Graph{Undirected}, states::Vector{Int}; t_var_name::String="t", p_var_name::String="p")
  @req length(states) == n_vertices(G) "need one and only one random variable in MarkovRing for each graph vertex"
  DiscreteGraphicalModel(G, states, Dict{Symbol, VarName}(:t => t_var_name, :p => p_var_name))
end

varnames(M::DiscreteGraphicalModel) = M.varnames
states(M::DiscreteGraphicalModel) = M.states
n_states(M::DiscreteGraphicalModel) = length(states(M))

function Base.show(io::IO, M::DiscreteGraphicalModel{T, L}) where {T, L}
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

@doc raw"""
    model_ring(M::DiscreteGraphicalModel; cached=false)

Return the ring in which the statistical model lives, together with a
Dict for indexing its generators. This is the same for directed and for
undirected graphs.

``` jldoctest 
julia> M = discrete_graphical_model(graph_from_edges([[1,2], [2,3]]), [2,2,2])
Discrete Graphical Model on a Undirected graph with 3 nodes and 2 edges with states [2, 2, 2]

julia> _, MR_gens = model_ring(M);

julia> MR_gens[1, 2, 1]
p[1,2,1]

julia> MR_gens[(1, 2, 1)]
p[1,2,1]
```
"""
@attr Tuple{
  ModelRing,
  GenDict
} function model_ring(M::DiscreteGraphicalModel; cached=false)
  random_variables = 1:n_states(M)
  state_spaces = [1:s for s in states(M)]
  varindices = collect(Iterators.product(state_spaces...))

  # TODO: the base ring should come from the model, leaving as QQ for now
  return model_ring(QQ, varnames(M)[:p] => varindices)
end

function state_space(M::DiscreteGraphicalModel, C::Vector{Int})
  return length(C) == 1 ? (1:states(M)[C[1]]) : Iterators.product([1:states(M)[i] for i in C]...)
end

state_space(M::DiscreteGraphicalModel, C::Set{Int}) = state_space(M, collect(C))

@doc raw"""
    parameter_ring(M::DiscreteGraphicalModel; cached=false)

Return the ring of parameters of the statistical model, together with a
Dict for indexing its generators.

``` jldoctest
julia> G = graph_from_edges([[1,2], [2,3]]);

julia> M = discrete_graphical_model(G, [2,2,2])
Discrete Graphical Model on a Undirected graph with 3 nodes and 2 edges with states [2, 2, 2]

julia> _, PR_gens = parameter_ring(M);

julia> PR_gens[[1, 2], (2, 2)]
t{1,2}(2, 2)

julia> C = maximal_cliques(G);

julia> PR_gens[first(C), (1, 2)]
t{2,3}(1, 2)
```
"""
@attr Tuple{
  MPolyRing,
  GenDict
} function parameter_ring(M::DiscreteGraphicalModel{Graph{Undirected}, T}; cached=false) where T
  cliques = maximal_cliques(graph(M))
  Xs = [state_space(M, C) for C in cliques]
  params = [(C, x) for (C, X) in Iterators.zip(cliques, Xs) for x in X]
  gen_names = [varnames(M)[:t] * "{" *join(string.(sort(collect(C))), ",") * "}" * string(x) for (C, x) in params]
  R, t = polynomial_ring(QQ, gen_names; cached=cached)
  gens_dict = Dict(zip(params, t))
  return (R, gens_dict)
end

@doc raw"""
    vanishing_ideal(M::DiscreteGraphicalModel{Graph{Undirected}, L} where L

The vanishing ideal is toric of an undirected discrete graphical model
is toric. Therefore its `vanishing_ideal` method defaults to `algorithm=:markov`.
"""
function vanishing_ideal(M::DiscreteGraphicalModel{Graph{Undirected}, L}) where L
  vanishing_ideal(M; algorithm=:markov)
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
Discrete Graphical Model on a Undirected graph with 3 nodes and 2 edges with states [2, 2, 2]

julia> parametrization(M)
Ring homomorphism
  from multivariate polynomial ring in 8 variables over QQ
  to multivariate polynomial ring in 8 variables over QQ
defined by
  p[1,1,1] -> t[1, 2](1, 1)*t[2, 3](1, 1)
  p[2,1,1] -> t[1, 2](2, 1)*t[2, 3](1, 1)
  p[1,2,1] -> t[1, 2](1, 2)*t[2, 3](2, 1)
  p[2,2,1] -> t[1, 2](2, 2)*t[2, 3](2, 1)
  p[1,1,2] -> t[1, 2](1, 1)*t[2, 3](1, 2)
  p[2,1,2] -> t[1, 2](2, 1)*t[2, 3](1, 2)
  p[1,2,2] -> t[1, 2](1, 2)*t[2, 3](2, 2)
  p[2,2,2] -> t[1, 2](2, 2)*t[2, 3](2, 2)
```
"""
function parametrization(M::DiscreteGraphicalModel{Graph{Undirected}, L}) where L
  S, pd = model_ring(M)
  R, td = parameter_ring(M)
  images = elem_type(R)[]
  for p in gens(_ring(S))
    s = S[p] # get label of the variable
    k = first(keys(td))
    push!(images, prod(td[k] for k in keys(td) if k[2] == s[collect(k[1])]))
  end
  hom(S, R, images)
end

@doc raw"""
    ci_structure(M::DiscreteGraphicalModel{Graph{Undirected}, L})

Return the set of elementary CI statements which are satisfied by every distribution in the
graphical model. This is the same as the `global_markov` property of the underlying graph.
"""
function ci_structure(M::DiscreteGraphicalModel{Graph{Undirected}, L}) where L
  global_markov(graph(M))
end

###################################################################################
#
#       Directed Discrete Graphical Models
#
###################################################################################

@doc raw"""
    discrete_graphical_model(G::Graph{Directed}, states::Vector{Int}; q_var_name::String="q")

The parametric statistical model associated to a directed acyclic graph.
It contains a directed acylic graph `G`, a MarkovRing `S` where the vanishing ideal of the model naturally lives,
and a parameter ring whose variables `q[i][j](j_1, j_2, ..., j_#pa(i))` correspond to the conditional distribution
of the node i given its parents.

``` jldoctest 
julia> M = discrete_graphical_model(graph_from_edges(Directed, [[1,3], [2,3], [3,4]]), [2,2,2,2])
Discrete Graphical Model on a Directed graph with 4 nodes and 3 edges with states [2, 2, 2, 2]
```
"""
function discrete_graphical_model(G::Graph{Directed}, states::Vector{Int}; q_var_name::String="q", p_var_name::String="p")
  @req length(states) == n_vertices(G) "need one and only one random variable in MarkovRing for each graph vertex"
  @req is_acyclic(G) "$G must be an acyclic graph"
  DiscreteGraphicalModel(G, states, Dict{Symbol, VarName}(:q => q_var_name, :p => p_var_name))
end

@doc raw"""
    parameter_ring(M::DiscreteGraphicalModel{Graph{Directed}, T}; cached=false)

Return the ring of parameters of the statistical model, together with a
Dict for indexing its generators.

There is a "hidden" homogenizing variable named "_h" in the ring. It can be
obtained via `last(gens(R))` but is not listed in the dictionary of generators.

``` jldoctest 
julia> M = discrete_graphical_model(graph_from_edges(Directed, [[1,3], [2,3], [3,4]]), [2,2,2,2])
Discrete Graphical Model on a Directed graph with 4 nodes and 3 edges with states [2, 2, 2, 2]

julia> _, PR_gens = parameter_ring(M);

julia> PR_gens[3, 2, [2]]
q[3](2 | 2)

julia> PR_gens[1, 1, Set{Int}()]
q[1](1)
```
"""
@attr Tuple{
  MPolyRing,
  GenDict
} function parameter_ring(M::DiscreteGraphicalModel{Graph{Directed}, T}; cached=false) where T
  G = graph(M)
  params = Tuple{Int, Int, Set{Int}}[]
  for i in 1:n_vertices(G)
    pa = parents(G, i);
    if length(pa) == 0
      append!(params, (i, x, Set{Int}([])) for x in state_space(M, [i]))
    else
      append!(params, (i, x, Set{Int}(Iterators.flatten([y]))) for x in state_space(M, [i]) for y in state_space(M, pa))
    end
  end
  gen_names = [varnames(M)[:q] * "[" * string(i) * "](" * string(x) * (length(y) > 0 ? " | " * join(y, ", ") : "") * ")" for (i, x, y) in params]
  push!(gen_names, "_h")
  R, q = polynomial_ring(QQ, gen_names; cached=cached);
  gens_dict = Dict(p => q[k] for (k, p) in enumerate(params))
  return (R, gens_dict)
end

@doc raw"""
  parameterization(M::GraphicalModel{Graph{Directed}, MarkovRing})

Creates the polynomial map which parameterizes the vanishing ideal of the directed discrete graphical model `M`.
Every probability coordinate in the `model_ring` is expressed as the product of conditional densities of each
node given its parents. The conditional probability parameters are subject to linear constraints.

## Examples

``` jldoctest
julia> M = discrete_graphical_model(graph_from_edges(Directed, [[1,3], [2,3]]), [2,2,2])
Discrete Graphical Model on a Directed graph with 3 nodes and 2 edges with states [2, 2, 2]

julia> phi = parametrization(M)
Ring homomorphism
  from multivariate polynomial ring in 8 variables over QQ
  to multivariate polynomial ring in 13 variables over QQ
defined by
  p[1,1,1] -> q[1](1)*q[2](1)*_h - q[1](1)*q[3](2 | 1)*_h - q[1](2)*q[2](2)*q[3](2 | 1) + q[2](2)*q[3](2 | 1)*_h
  p[2,1,1] -> q[1](2)*q[2](1)*_h + q[1](2)*q[2](2)*q[3](2 | 2, 1) - q[1](2)*q[3](2 | 2, 1)*_h
  p[1,2,1] -> q[1](1)*q[2](2)*_h + q[1](2)*q[2](2)*q[3](2 | 2, 1) - q[2](2)*q[3](2 | 2, 1)*_h
  p[2,2,1] -> -q[1](2)*q[2](2)*q[3](2 | 2) + q[1](2)*q[2](2)*_h
  p[1,1,2] -> q[1](1)*q[3](2 | 1)*_h + q[1](2)*q[2](2)*q[3](2 | 1) - q[2](2)*q[3](2 | 1)*_h
  p[2,1,2] -> -q[1](2)*q[2](2)*q[3](2 | 2, 1) + q[1](2)*q[3](2 | 2, 1)*_h
  p[1,2,2] -> -q[1](2)*q[2](2)*q[3](2 | 2, 1) + q[2](2)*q[3](2 | 2, 1)*_h
  p[2,2,2] -> q[1](2)*q[2](2)*q[3](2 | 2)

julia> I = kernel(phi)
Ideal generated by
  -p[1,1,1]*p[2,2,1] - p[1,1,1]*p[2,2,2] + p[2,1,1]*p[1,2,1] + p[2,1,1]*p[1,2,2] + p[1,2,1]*p[2,1,2] - p[2,2,1]*p[1,1,2] - p[1,1,2]*p[2,2,2] + p[2,1,2]*p[1,2,2]
```
"""
function parametrization(M::DiscreteGraphicalModel{Graph{Directed}, T}) where T
  G = graph(M)
  n = n_vertices(G)
  S, pd = model_ring(M)
  R, qd = parameter_ring(M)
  h = last(gens(R))
  images = elem_type(R)[]
  for p in gens(_ring(S))
    s = collect(S[p]) # get label of the variable but as a vector for type reasons
    push!(images, prod(qd[i, s[i], s[parents(G, i)]] for i in 1:n))
  end
  # There are linear constraints to enforce on the `images`. Instead of
  # doing this by hand, we employ a Gröbner basis of the linear ideal.
  condprob = elem_type(S)[]
  for i in 1:n
    pa = parents(G, i)
    if length(pa) == 0
      push!(condprob, h - sum(qd[i, x, Int64[]] for x in state_space(M, [i])))
    else
      for y in state_space(M, pa)
        push!(condprob, h - sum(qd[i, x, collect(Iterators.flatten([y]))] for x in state_space(M, [i])))
      end
    end
  end
  I = ideal(condprob)
  return hom(S, R, map(Base.Fix2(reduce, gens(I)), images))
end

@doc raw"""
    ci_structure(M::DiscreteGraphicalModel{Graph{Directed}, L})

Return the set of elementary CI statements which are satisfied by every distribution in the
graphical model. This is the same as the `global_markov` property of the underlying graph.
"""
function ci_structure(M::DiscreteGraphicalModel{Graph{Directed}, L}) where L
  global_markov(graph(M))
end
