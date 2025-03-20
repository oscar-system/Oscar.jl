###################################################################################
#
#       (Un)Directed Gaussian Graphical Models
#
###################################################################################

@attributes mutable struct GaussianGraphicalModel{T, L} <: GraphicalModel{T, L}
  graph::Graph{T}
  labelings::L
  varnames::Vector{VarName}
  function GaussianGraphicalModel(G::Graph{T}, var_names::Vector{VarName}) where T <: GraphTypes
    graph_maps = NamedTuple(_graph_maps(G))
    graph_maps = isempty(graph_maps) ? nothing : graph_maps
    return new{T, typeof(graph_maps)}(G, graph_maps, var_names)
  end
end

@doc raw"""
     gaussian_graphical_model(G::Graph{Directed}; s_varname::VarName="s", l_varname::VarName="l", w_varname::VarName="w")
     gaussian_graphical_model(G::Graph{Undirected}; s_varname::VarName="s", k_varname::VarName="k")

A parametric statistical model associated to a graph `G`.
If `cached` is `true`, the internally generated polynomial ring will be cached.

## Examples

```jldoctest
julia> M = gaussian_graphical_model(graph_from_edges(Directed, [[1,2], [2,3]]))
Gaussian Graphical Model on a Directed graph with 3 nodes and 2 edges

julia> graph(M)
Directed graph with 3 nodes and the following edges:
(1, 2)(2, 3)

```
"""
function gaussian_graphical_model(G::Graph{Directed}; s_varname::VarName="s", l_varname::VarName="l", w_varname::VarName="w")
  GaussianGraphicalModel(G, VarName[s_varname, l_varname, w_varname])
end

function gaussian_graphical_model(G::Graph{Undirected}; s_varname::VarName="s", k_varname::VarName="k")
  GaussianGraphicalModel(G, VarName[s_varname, k_varname])
end

function Base.show(io::IO, M::GaussianGraphicalModel{T, L}) where {T, L}
  io = pretty(io)
  if is_terse(io)
    print(io, "Gaussian Graphical Model on a")
    !(L == Nothing) && print(io, " labelled")
    print(io, " $T graph")
  else
    print(io, "Gaussian Graphical Model on a $(graph(M))")
  end
end

@attr function probability_ring(GM::GaussianGraphicalModel; cached=false)
  n = n_vertices(graph(GM))
  varindices = [(i, j) for i in 1:n for j in i:n]
  varnames = ["$(GM.varnames[1])[$(i), $(j)]" for (i, j) in varindices]
  polynomial_ring(QQ, varnames; cached=false)[1]
end

@doc raw"""
    covariance_matrix(GM::GaussianGraphicalModel)

Return the covariance matrix associated to the graphical model `GM` as a matrix over the underlying probability ring of `GM`,
see [`probability_ring`](@ref)

## Examples

```jldoctest
julia> GM = gaussian_graphical_model(graph_from_edges(Directed, [[1,2], [2,3]]))
Gaussian Graphical Model on a Directed graph with 3 nodes and 2 edges

julia> covariance_matrix(GM)
[s[1, 1]   s[1, 2]   s[1, 3]]
[s[1, 2]   s[2, 2]   s[2, 3]]
[s[1, 3]   s[2, 3]   s[3, 3]]

```
"""
function covariance_matrix(GM::GaussianGraphicalModel)
  S = probability_ring(GM)
  cov_mat = upper_triangular_matrix(gens(S))
  n = n_vertices(graph(GM))
  # turn upper triangular matrix into a symmetric one
  for i in 1:n
    for j in i + 1:n
      cov_mat[j, i] = cov_mat[i, j]
    end
  end
  return cov_mat
end

###################################################################################
#
#       Directed parametrization
#
###################################################################################

@attr function parameter_ring(GM::GaussianGraphicalModel{Directed, T}; cached=false) where T
  G = graph(GM)
  varnames = (["$(GM.varnames[2])[$(src(e)), $(dst(e))]" for e in edges(G)], ["$(GM.varnames[3])[$(v)]" for v in vertices(G)])
  polynomial_ring(QQ, varnames; cached=cached)[1]
end

function parameter_ring_gens(GM::GaussianGraphicalModel{Directed, T}) where T
  G = graph(GM)
  R = parameter_ring(GM)
  varnames = (["$(GM.varnames[2])[$(src(e)), $(dst(e))]" for e in edges(G)],
              ["$(GM.varnames[3])[$(v)]" for v in vertices(G)])
  edge_vars, vertex_vars = AbstractAlgebra.reshape_to_varnames(gens(R), varnames)
  
  merge(Dict(e => edge_vars[i] for (i, e) in enumerate(edges(G))),
        Dict(v => vertex_vars[v] for v in 1:n_vertices(G)))
end


@doc raw"""
    directed_edges_matrix(M::GaussianGraphicalModel{Graph{Directed}, T}) where T

Create the weighted adjacency matrix $\Lambda$ of a directed graph `G` whose entries are the parameter ring of the graphical model `M`.

## Examples

```jldoctest
julia> M = graphical_model(graph_from_edges(Directed, [[1,2], [2,3]]), gaussian_ring(3))
Gaussian graphical model on a directed graph with edges:
(1, 2), (2, 3)

julia> directed_edges_matrix(M)
[0   l[1, 2]         0]
[0         0   l[2, 3]]
[0         0         0]
```
"""
function directed_edges_matrix(M::GaussianGraphicalModel{Directed, T}) where T
  G = graph(M)
  gens_dict = parameter_ring_gens(M)
  n =  n_vertices(G)
  L = zero_matrix(parameter_ring(M), n, n)
  for e in edges(G)
    L[src(e), dst(e)] = gens_dict[e]
  end
  return L
end

@doc raw"""
    error_covariance_matrix(M::GraphicalModel{Graph{Directed}, T}) where T

Create the covariance matrix $ \Omega $ of the independent error terms in a directed Gaussian graphical model `M`

## Examples

```jldoctest
julia> M = graphical_model(graph_from_edges(Directed, [[1,2], [2,3]]), gaussian_ring(3))
Gaussian graphical model on a directed graph with edges:
(1, 2), (2, 3)

julia> error_covariance_matrix(M)
[w[1]      0      0]
[   0   w[2]      0]
[   0      0   w[3]]

```
"""
function error_covariance_matrix(M::GaussianGraphicalModel{Directed, T}) where T
  G = graph(M)
  gens_dict = parameter_ring_gens(M)
  W = diagonal_matrix(parameter_ring(M), [gens_dict[i] for i in 1:n_vertices(G)])
end

@doc raw"""
    parametrization(M::GaussianGraphicalModel{Directed, T}) where T

Create the polynomial map which parametrizes the vanishing ideal of the directed Gaussian graphical model `M`.
The vanishing ideal of the statistical model is the kernel of this map. This ring map is the pull back of the parametrization $\phi_G$ given by
$(Id - \Lambda)^{-T} \Omega (Id - \Lambda)^{T} \mapsto \Sigma$ where $\Lambda =$ `directed_edges_matrix(M)` and  $ \Omega = $  `error_covariance_matrix(M)`.

## Examples

```jldoctest
julia> M = graphical_model(graph_from_edges(Directed, [[1,2], [2,3]]), gaussian_ring(3))
Gaussian graphical model on a directed graph with edges:
(1, 2), (2, 3)

julia> parametrization(M)
Ring homomorphism
  from multivariate polynomial ring in 6 variables over QQ
  to multivariate polynomial ring in 5 variables over QQ
defined by
  s[1, 1] -> w[1]
  s[1, 2] -> l[1, 2]*w[1]
  s[1, 3] -> l[1, 2]*l[2, 3]*w[1]
  s[2, 2] -> l[1, 2]^2*w[1] + w[2]
  s[2, 3] -> l[1, 2]^2*l[2, 3]*w[1] + l[2, 3]*w[2]
  s[3, 3] -> l[1, 2]^2*l[2, 3]^2*w[1] + l[2, 3]^2*w[2] + w[3]
```
"""
function parametrization(M::GaussianGraphicalModel{Directed, T}) where T
  S = probability_ring(M)
  R = parameter_ring(M)
  G = graph(M)
  L = directed_edges_matrix(M)
  Id = identity_matrix(R, n_vertices(G))
  W = error_covariance_matrix(M)
  Sigma = transpose(inv(Id - L)) * W * inv(Id - L)
  hom(S, R, reduce(vcat, [[Sigma[i,j] for j in i:n_vertices(G)] for i in 1:n_vertices(G)]))
end

###################################################################################
#
#       Undirected parametrization
#
###################################################################################
@attr function parameter_ring(GM::GaussianGraphicalModel{Undirected, T}; cached=false) where T
  G = graph(GM)
  varnames = (["$(GM.varnames[2])[$(src(e)), $(dst(e))]" for e in edges(G)], ["$(GM.varnames[2])[$(v), $(v)]" for v in vertices(G)])
  polynomial_ring(QQ, varnames; cached=cached)[1]
end

function parameter_ring_gens(GM::GaussianGraphicalModel{Undirected, T}) where T
  G = graph(GM)
  R = parameter_ring(GM)
  varnames = (["$(GM.varnames[2])[$(src(e)), $(dst(e))]" for e in edges(G)],
              ["$(GM.varnames[2])[$(v), $(v)]" for v in vertices(G)])
  edge_vars, vertex_vars = AbstractAlgebra.reshape_to_varnames(gens(R), varnames)
  
  merge(Dict(e => edge_vars[i] for (i, e) in enumerate(edges(G))),
        Dict(v => vertex_vars[v] for v in 1:n_vertices(G)))
end

@doc raw"""
    concentration_matrix(M::GaussianGraphicalModel{Undirected, T}) where T

Create the concentration matrix `K` of an undirected Gaussian graphical model which is a symmetric positive definite matrix
whose nonzero entries correspond to the edges of the associated graph.

## Examples

```jldoctest
julia> M = gaussian_graphical_model(graph_from_edges([[1,2], [2,3]]))
Gaussian graphical model on an undirected graph with edges:
(1, 2), (2, 3)

julia> concentration_matrix(M)
[k[1, 1]   k[1, 2]         0]
[k[1, 2]   k[2, 2]   k[2, 3]]
[      0   k[2, 3]   k[3, 3]]
```
"""
function concentration_matrix(M::GaussianGraphicalModel{Undirected, T}) where T
  G = graph(M)
  k = parameter_ring_gens(M)
  K = diagonal_matrix(parameter_ring(M), [k[i] for i in 1:n_vertices(G)])

  for e in edges(G)
    K[dst(e), src(e)] = k[e]
    K[src(e), dst(e)] = k[e]
  end
  return K
end

@doc raw"""
    parametrization(M::GaussianGraphicalModel{Undirected, Nothing})

Create the polynomial map which parametrizes the vanishing ideal of the undirected Gaussian graphical model `M`.
The vanishing ideal of the statistical model is the kernel of this map. This ring map is the pull back of the parametrization $\phi_G$ given by
$ K \mapsto K^{-1}$ where $ K = $  `concentration_matrix(M)` and the entries of $ K^{-1} $ are given by the standard cofactor formula.

## Examples

```jldoctest
julia> M = graphical_model(graph_from_edges([[1,2], [2,3]]), gaussian_ring(3))
Gaussian graphical model on an undirected graph with edges:
(1, 2), (2, 3)

julia> parametrization(M)
Ring homomorphism
  from multivariate polynomial ring in 6 variables over QQ
  to fraction field of multivariate polynomial ring
defined by
  s[1, 1] -> (k[2, 2]*k[3, 3] - k[2, 3]^2)//(k[1, 1]*k[2, 2]*k[3, 3] - k[1, 1]*k[2, 3]^2 - k[1, 2]^2*k[3, 3])
  s[1, 2] -> (-k[1, 2]*k[3, 3])//(k[1, 1]*k[2, 2]*k[3, 3] - k[1, 1]*k[2, 3]^2 - k[1, 2]^2*k[3, 3])
  s[1, 3] -> (k[1, 2]*k[2, 3])//(k[1, 1]*k[2, 2]*k[3, 3] - k[1, 1]*k[2, 3]^2 - k[1, 2]^2*k[3, 3])
  s[2, 2] -> (k[1, 1]*k[3, 3])//(k[1, 1]*k[2, 2]*k[3, 3] - k[1, 1]*k[2, 3]^2 - k[1, 2]^2*k[3, 3])
  s[2, 3] -> (-k[1, 1]*k[2, 3])//(k[1, 1]*k[2, 2]*k[3, 3] - k[1, 1]*k[2, 3]^2 - k[1, 2]^2*k[3, 3])
  s[3, 3] -> (k[1, 1]*k[2, 2] - k[1, 2]^2)//(k[1, 1]*k[2, 2]*k[3, 3] - k[1, 1]*k[2, 3]^2 - k[1, 2]^2*k[3, 3])
```
"""
function parametrization(M::GaussianGraphicalModel{Undirected, T}) where T
  G = graph(M)
  S = probability_ring(M)
  R = parameter_ring(M)
  K = concentration_matrix(M)
  adj = adjugate(K)
  hom(S, fraction_field(R), reduce(vcat, [[adj[i,j]//det(K) for j in i:n_vertices(G)] for i in 1:n_vertices(G)]))
end

@doc raw"""
    vanishing_ideal(M::GaussianGraphicalModel{Undirected, Nothing})

Compute the vanishing ideal of the undirected Gaussian graphical model `M`.
This is done by saturating the ideal given by $ \Sigma K - Id $ by the determinant of $ K $
and then eliminating all variables `k[i,j]` where $ K =$ `concentration_matrix(M)`.

## Examples

```jldoctest
julia> M = gaussian_graphical_model(graph_from_edges([[1,2], [2,3]]))
Gaussian Graphical Model on a Undirected graph with 3 nodes and 2 edges

julia> vanishing_ideal(M)
Ideal generated by
  -s[1, 2]*s[2, 3] + s[1, 3]*s[2, 2]
```
"""
function vanishing_ideal(M::GaussianGraphicalModel{Undirected, T}) where T
  G = graph(M)
  S = probability_ring(M)
  R = parameter_ring(M)

  # there simply must be a better way to do this but I cannot find one currently
  elim_ring, elim_gens = polynomial_ring(coefficient_ring(R), vcat([string(x) for x in gens(R)], [string(x) for x in gens(S)]))
  inject_R = hom(R, elim_ring, elim_gens[1:ngens(R)])
  inject_S = hom(S, elim_ring, elim_gens[ngens(R)+1:ngens(elim_ring)])

  Sigma = inject_S.(covariance_matrix(M))
  K = inject_R.(concentration_matrix(M))

  elim_ideal = ideal(reduce(vcat, Sigma*K - identity_matrix(elim_ring, n_vertices(G))))
  eliminate(saturation(elim_ideal, ideal([det(K)])), elim_gens[1:ngens(R)])
end
