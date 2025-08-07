# this Any here makes me uncomfortable
const GenDictType = Dict{Any, QQMPolyRingElem}
###################################################################################
#
#       (Un)Directed Gaussian Graphical Models
#
###################################################################################

@attributes mutable struct GaussianGraphicalModel{T, L} <: GraphicalModel{T, L}
  graph::Graph{T}
  labelings::L
  varnames::Dict{Symbol,VarName}
  function GaussianGraphicalModel(G::Graph{T},
                                  var_names::Dict{Symbol, VarName}) where T <: GraphTypes
    graph_maps = NamedTuple(_graph_maps(G))
    graph_maps = isempty(graph_maps) ? nothing : graph_maps
    return new{T, typeof(graph_maps)}(G, graph_maps, var_names)
  end

  function GaussianGraphicalModel(G::MixedGraph, var_names::Vector{VarName})
    #TODO figure out how to deal with labelings on MixedGraphs
    # for now just use Nothing
    return new{Mixed, Nothing}(G, nothing, var_names)
  end
end

@doc raw"""
    gaussian_graphical_model(G::Graph{Directed}; s_varname::VarName=:s, l_varname::VarName=:ll, w_varname::VarName=:w)
    gaussian_graphical_model(G::Graph{Undirected}; s_varname::VarName=:s, k_varname::VarName=:k)
Given a graph `G` construct a gaussian graphical model for `G`. Optionally one can set the letter used for the variable of the covariance matrix by setting `s_varname`. 

## Examples

```jldoctest
julia> G = graph_from_edges(Directed, [[1, 2], [2, 3]])
Directed graph with 3 nodes and the following edges:
(1, 2)(2, 3)

julia> GM = gaussian_graphical_model(G)
Gaussian Graphical Model on a Directed graph with 3 nodes and 2 edges

julia> typeof(GM)
GaussianGraphicalModel{Directed, Nothing}

julia> CG = graph_from_labeled_edges(Dict((1, 2) => "red", (2, 3) => "green"); name=:color)
Undirected graph with 3 nodes and the following labeling(s):
label: color
(2, 1) -> red
(3, 2) -> green

julia> GM = gaussian_graphical_model(CG)
Gaussian Graphical Model on a Undirected graph with 3 nodes and 2 edges with labeling(s) [:color]

julia> typeof(GM)
GaussianGraphicalModel{Undirected, @NamedTuple{color::Oscar.GraphMap{Undirected, Polymake.LibPolymake.EdgeMapAllocated{Undirected, CxxWrap.StdLib.StdString}, Nothing}}}
```
"""
function gaussian_graphical_model(G::Graph{Directed}; s_varname::VarName="s", l_varname::VarName="l", w_varname::VarName="w")
  GaussianGraphicalModel(G, Dict{Symbol, VarName}(:s => s_varname, :l => l_varname, :w => w_varname))
end

function gaussian_graphical_model(G::Graph{Undirected}; s_varname::VarName="s", k_varname::VarName="k")
  GaussianGraphicalModel(G, Dict{Symbol, VarName}(:s => s_varname, :k => k_varname))
end

varnames(GM::GaussianGraphicalModel) = GM.varnames

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

#TODO do the gens returned here need to be a dict? we'll need to be consisted across all models
@attr Tuple{
  MPolyRing,
  Vector{QQMPolyRingElem}
} function model_ring(GM::GaussianGraphicalModel; cached=false)
  n = n_vertices(graph(GM))
  varindices = [(i, j) for i in 1:n for j in i:n]
  gen_names = ["$(GM.varnames[:s])[$(i), $(j)]" for (i, j) in varindices]
  # we'll need to align on the return value of the gens, it might need to be Dict
  return polynomial_ring(QQ, gen_names; cached=false)
end

@doc raw"""
    covariance_matrix(GM::GaussianGraphicalModel)

Return the covariance matrix associated to the graphical model `GM` as a matrix over the underlying model ring of `GM`,
see [`model_ring`](@ref)

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
  _, s = model_ring(GM)
  cov_mat = upper_triangular_matrix(s)
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

@attr Tuple{
  MPolyRing,
  GenDictType
} function parameter_ring(GM::GaussianGraphicalModel{Directed, T}; cached=false) where T
  G = graph(GM)
  gen_names = (["$(varnames(GM)[:l])[$(src(e)), $(dst(e))]" for e in edges(G)],
              ["$(varnames(GM)[:w])[$(v)]" for v in vertices(G)])
  R, e_gens, v_gens = polynomial_ring(QQ, gen_names; cached=cached)
  gens_dict = merge(Dict(e => e_gens[i] for (i, e) in enumerate(edges(G))),
                    Dict(v => v_gens[v] for v in 1:n_vertices(G)))
  return R, gens_dict
end

@doc raw"""
    directed_edges_matrix(M::GaussianGraphicalModel{Directed, L}) where L

Create the weighted adjacency matrix $\Lambda$ of a directed graph `G` whose entries are the parameter ring of the graphical model `M`.

## Examples

```jldoctest
julia> GM = gaussian_graphical_model(graph_from_edges(Directed, [[1,2], [2,3]]))
Gaussian Graphical Model on a Directed graph with 3 nodes and 2 edges

julia> directed_edges_matrix(GM)
[0   l[1, 2]         0]
[0         0   l[2, 3]]
[0         0         0]

```
"""
function directed_edges_matrix(M::GaussianGraphicalModel{Directed, L}) where L
  G = graph(M)
  R, gens_dict = parameter_ring(M)
  n =  n_vertices(G)
  lambda = zero_matrix(R, n, n)
  for e in edges(G)
    lambda[src(e), dst(e)] = gens_dict[e]
  end
  return lambda
end

@doc raw"""
    error_covariance_matrix(M::GaussianGraphicalModel{Directed, L}) where L

Create the covariance matrix $ \Omega $ of the independent error terms in a directed Gaussian graphical model `M`

## Examples

```jldoctest
julia> M = gaussian_graphical_model(graph_from_edges(Directed, [[1,2], [2,3]]))
Gaussian Graphical Model on a Directed graph with 3 nodes and 2 edges

julia> error_covariance_matrix(M)
[w[1]      0      0]
[   0   w[2]      0]
[   0      0   w[3]]

```
"""
function error_covariance_matrix(M::GaussianGraphicalModel{Directed, L}) where L
  G = graph(M)
  R, gens_dict = parameter_ring(M)
  W = diagonal_matrix(R, [gens_dict[i] for i in 1:n_vertices(G)])
end

@doc raw"""
    parametrization(M::GaussianGraphicalModel{Directed, L}) where L

Create the polynomial map which parametrizes the vanishing ideal of the directed Gaussian graphical model `M`.
The vanishing ideal of the statistical model is the kernel of this map. This ring map is the pull back of the parametrization $\phi_G$ given by
$(Id - \Lambda)^{-T} \Omega (Id - \Lambda)^{T} \mapsto \Sigma$ where $\Lambda =$ `directed_edges_matrix(M)` and  $ \Omega = $  `error_covariance_matrix(M)`.

## Examples

```jldoctest
julia> M = gaussian_graphical_model(graph_from_edges(Directed, [[1,2], [2,3]]))
Gaussian Graphical Model on a Directed graph with 3 nodes and 2 edges

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
function parametrization(M::GaussianGraphicalModel{Directed, L}) where L
  S, _ = model_ring(M)
  R, _ = parameter_ring(M)
  G = graph(M)
  lambda = directed_edges_matrix(M)
  Id = identity_matrix(R, n_vertices(G))
  W = error_covariance_matrix(M)
  Sigma = transpose(inv(Id - lambda)) * W * inv(Id - lambda)
  hom(S, R, reduce(vcat, [[Sigma[i,j] for j in i:n_vertices(G)] for i in 1:n_vertices(G)]))
end

###################################################################################
#
#       Undirected parametrization
#
###################################################################################
@attr Tuple{MPolyRing, GenDictType} function parameter_ring(GM::GaussianGraphicalModel{Undirected, T}; cached=false) where T
  G = graph(GM)
  gen_names = (["$(varnames(GM)[:k])[$(src(e)), $(dst(e))]" for e in edges(G)],
              ["$(varnames(GM)[:k])[$(v), $(v)]" for v in vertices(G)])
  R, (e_gens, v_gens) = polynomial_ring(QQ, gen_names; cached=cached)[1]
  gens_dict = merge(Dict(e => edge_vars[i] for (i, e) in enumerate(edges(G))),
                    Dict(v => vertex_vars[v] for v in 1:n_vertices(G)))
  return R, gens_dict
end

@doc raw"""
    concentration_matrix(M::GaussianGraphicalModel{Undirected, T}) where T

Create the concentration matrix `K` of an undirected Gaussian graphical model which is a symmetric positive definite matrix
whose nonzero entries correspond to the edges of the associated graph.

## Examples

```jldoctest
julia> M = gaussian_graphical_model(graph_from_edges([[1,2], [2,3]]))
Gaussian Graphical Model on a Undirected graph with 3 nodes and 2 edges

julia> concentration_matrix(M)
[k[1, 1]   k[2, 1]         0]
[k[2, 1]   k[2, 2]   k[3, 2]]
[      0   k[3, 2]   k[3, 3]]

```
"""
function concentration_matrix(M::GaussianGraphicalModel{Undirected, T}) where T
  G = graph(M)
  R, k = parameter_ring(M)
  K = diagonal_matrix(R, [k[i] for i in 1:n_vertices(G)])

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
julia> M = gaussian_graphical_model(graph_from_edges([[1,2],[2,3]]))
Gaussian Graphical Model on a Undirected graph with 3 nodes and 2 edges

julia> parametrization(M)
Ring homomorphism
  from multivariate polynomial ring in 6 variables over QQ
  to localization of multivariate polynomial ring in 5 variables over QQ at products of (-k[2, 1]^2*k[3, 3] - k[3, 2]^2*k[1, 1] + k[1, 1]*k[2, 2]*k[3, 3])
defined by
  s[1, 1] -> (k[3, 2]^2 - k[2, 2]*k[3, 3])/(k[2, 1]^2*k[3, 3] + k[3, 2]^2*k[1, 1] - k[1, 1]*k[2, 2]*k[3, 3])
  s[1, 2] -> k[2, 1]*k[3, 3]/(k[2, 1]^2*k[3, 3] + k[3, 2]^2*k[1, 1] - k[1, 1]*k[2, 2]*k[3, 3])
  s[1, 3] -> -k[2, 1]*k[3, 2]/(k[2, 1]^2*k[3, 3] + k[3, 2]^2*k[1, 1] - k[1, 1]*k[2, 2]*k[3, 3])
  s[2, 2] -> -k[1, 1]*k[3, 3]/(k[2, 1]^2*k[3, 3] + k[3, 2]^2*k[1, 1] - k[1, 1]*k[2, 2]*k[3, 3])
  s[2, 3] -> k[3, 2]*k[1, 1]/(k[2, 1]^2*k[3, 3] + k[3, 2]^2*k[1, 1] - k[1, 1]*k[2, 2]*k[3, 3])
  s[3, 3] -> (k[2, 1]^2 - k[1, 1]*k[2, 2])/(k[2, 1]^2*k[3, 3] + k[3, 2]^2*k[1, 1] - k[1, 1]*k[2, 2]*k[3, 3])
```
"""
function parametrization(M::GaussianGraphicalModel{Undirected, T}) where T
  G = graph(M)
  S = model_ring(M)[1]
  R = parameter_ring(M)[1]
  K = concentration_matrix(M)
  Rloc, iota = localization(R, powers_of_element(det(K)))
  adj = adjugate(K)
  hom(S, Rloc, reduce(vcat, [[iota(adj[i,j])/iota(det(K)) for j in i:n_vertices(G)] for i in 1:n_vertices(G)]))
end

###################################################################################
#
#      Gaussian Rings
#
###################################################################################

struct GaussianRing
  ring::MPolyRing
  #gens::Dict i removed this since all the info can be retrieved from the covariance matrix
  covariance_matrix::MatElem
end

@doc raw"""
    gaussian_ring(F::Field, n::Int; s_var_name::VarName=:s, cached=false)
    gaussian_ring(n::Int; s_var_name::VarName=:s, cached=false)
    gaussian_ring(GM::GaussianGraphicalmodel)

A polynomial ring whose variables correspond to the entries of a covariance matrix of `n` Gaussian random variables.
It is a multivariate polynomial ring whose variables are named `s[i,j]`and whose coefficient field is `QQ` by default or one can specify the field `F` by passing it as the first argument.


If `cached` is `true`, the internally generated polynomial ring will be cached.
## Examples
```jldoctest
julia> R = gaussian_ring(3)
Gaussian ring over Rational field in 6 variables
s[1, 1], s[1, 2], s[1, 3], s[2, 2], s[2, 3], s[3, 3]

julia> M = gaussian_graphical_model(graph_from_edges([[1,2], [2,3]]))
Gaussian Graphical Model on a Undirected graph with 3 nodes and 2 edges

julia> gaussian_ring(M)
Gaussian ring over Rational field in 6 variables
s[1, 1], s[1, 2], s[1, 3], s[2, 2], s[2, 3], s[3, 3]

```
"""
function gaussian_ring(F::Field, n::Int; s_var_name::VarName=:s, cached=false)
  varindices = [Tuple([i,j]) for i in 1:n for j in i:n]
  gen_names = ["$(s_var_name)[$(i), $(j)]" for (i,j) in varindices]
  S, s = polynomial_ring(F, gen_names; cached=cached)
  d = Dict([varindices[i] => s[i] for i in 1:length(varindices)])
  cov_matrix = matrix([[i < j ? d[i,j] : d[j,i] for j in 1:n] for i in 1:n])

  #TODO this should probably also return a pair
  GaussianRing(S, cov_matrix)
end

gaussian_ring(n::Int; s_var_name::VarName=:s, cached=false) = gaussian_ring(QQ, n; s_var_name=s_var_name, cached=cached)

gaussian_ring(GM::GaussianGraphicalModel) = GaussianRing(model_ring(GM), covariance_matrix(GM))

covariance_matrix(GR::GaussianRing) = GR.covariance_matrix

function Base.show(io::IO, R::GaussianRing)
  coeffs = base_ring(R.ring)
  k = ngens(R.ring)
  print(io, "Gaussian ring over $(coeffs) in $(k) variables", "\n", chop(string(gens(R.ring)), head = 16, tail = 1))
end

