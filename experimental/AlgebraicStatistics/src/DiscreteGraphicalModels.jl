###################################################################################
#
#       Undirected Discrete Graphical Models
#
###################################################################################


function bron_kerbosch(G::Graph{Undirected}, R, P, X)

  cliques = []

  if length(P) == 0 && length(X) == 0
      
      push!(cliques, R)
  end

  Q = P
  Y = X

  for v in Q

      new_R = vcat(R, [v])
      new_Q = [p for p in Q if p in neighbors(G, v)]
      new_Y = [p for p in Y if p in neighbors(G, v)]

      append!(cliques, bron_kerbosch(G, new_R, new_Q, new_Y))

      Q = filter(i -> i != v, P)
      Y = vcat(X, v)
  end

  return cliques
end

function maximal_cliques(G::Graph{Undirected})

  max_cliques = bron_kerbosch(G, Vector{Int}[], vertices(G), Vector{Int}[]);

  map(C -> Vector{Int64}(C), unique(sort.(max_cliques)))
end



@doc raw"""
  graphical_model(G::Graph{Undirected}, S::MarkovRing, k_var_name::String="k")

A parametric statistical model associated to an undirected graph.
It contains an undirected graph `G`, a MarkovRing `S` where the vanishing ideal of the model naturally lives, 
and a parameter ring whose variables `t[C](i_1, i_2, ..., i_#C)` correspond to potential functions of each clique. 

## Examples

``` jldoctest 
julia> M = graphical_model(graph_from_edges([[1,2], [2,3]]), markov_ring("A" => 1:2, "B" => 1:2, "X" => 1:2))
discrete graphical model on an undirected graph with edges:
(1, 2), (2, 3)
```
"""
function graphical_model(G::Graph{Undirected}, S::MarkovRing, t_var_name::String="t")

  cliques = maximal_cliques(G)
  rvs = random_variables(S)

  R, t = polynomial_ring(QQ, [t_var_name*string(C)*string(i) for C in cliques for i in collect(state_space(S, rvs[C]))])

  param_dict = Dict([(Tuple(C), Dict()) for C in cliques])

  i = 1
  
  for C in cliques

      for clique_state in collect(state_space(S, rvs[C]))
          
          param_dict[Tuple(C)][Tuple(clique_state)] = t[i]
          i += 1
      end
  end

  GraphicalModel(G, S, R, param_dict)
end


#TODO method for printing output of a discrete undirected graphical model


@doc raw"""
  parameterization(M::GraphicalModel{Graph{Undirected}, MarkovRing})

Creates the polynomial map which parameterizes the vanishing ideal of the undirected discrete graphical model `M`.   
The vanishing ideal of the statistical model is the kernel of this map. This ring map is the pull back of the parameterization $\phi_G$
where every joint probability is given by a product of the potential functions of the cliques. 
## Examples

``` jldoctest
julia> M = graphical_model(graph_from_edges([[1,2], [2,3]]), markov_ring(2,2,2))
discrete graphical model on an undirected graph with edges:
(1, 2), (2, 3)

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
function parametrization(M::GraphicalModel{Graph{Undirected}, MarkovRing})

  G = graph(M)
  cliques = maximal_cliques(G)
  S = ring(M)
  R = param_ring(M)
  t = param_gens(M)
 
  images = []

  for ind in state_space(S)

      push!(images, prod(map(C -> t[Tuple(C)][Tuple(ind[C])], cliques)))
  end

  hom(ring(S), R, images)
end


###################################################################################
#
#       Directed Discrete Graphical Models
#
###################################################################################


@doc raw"""
  graphical_model(G::Graph{Directed}, S::MarkovRing, q_var_name::String="q")

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
function graphical_model(G::Graph{Directed}, S::MarkovRing, q_var_name::String="q")

    rvs = random_variables(S)

    R, q = polynomial_ring(QQ, vcat([q_var_name*string([i, j])*string(par_state) for i in vertices(G) for j in state_space(S, rvs[i]) for par_state in parental_state_space(G, i, S)], ["h"]))

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
function parametrization(M::GraphicalModel{Graph{Directed}, MarkovRing})

  G = graph(M)
  V = vertices(G)
  S = ring(M)
  rvs = random_variables(S)
  R = param_ring(M)
  q = param_gens(M)

  images = []

  for ind in state_space(S)

      push!(images, prod(map(i-> q[i, ind[i]][Tuple(ind[inneighbors(G, i)])], V)))
  end

  lin_constraints  = []

  for i in V
    for par_state in parental_state_space(G, i, S)

      push!(lin_constraints, last(gens(R)) - sum(map(j -> q[i, j][Tuple(par_state)] , state_space(S, rvs[i]))))
    end
  end
  I = ideal(lin_constraints)

  return hom(ring(S), R, map(f -> reduce(f, gens(I)), images))
end
