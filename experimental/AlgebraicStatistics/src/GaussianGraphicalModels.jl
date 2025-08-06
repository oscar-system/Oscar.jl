###################################################################################
#
#       Gaussian Rings
#
###################################################################################

export GaussianRing
export gaussian_ring, ring, gens, covariance_matrix
export graphical_model, parametrization
export directed_edges_matrix, error_covariance_matrix
export concentration_matrix

struct GaussianRing
  ring::MPolyRing
  gens::Dict{Tuple{Int64, Int64}, <:MPolyRingElem}
  covariance_matrix::MatElem
end

@doc raw"""
    gaussian_ring(n::Int; s_var_name::VarName="s", K::Field=QQ, cached=false)

A polynomial ring whose variables correspond to the entries of a covariance matrix of `n` Gaussian random variables.
It is a multivariate polynomial ring whose variables are named `s[i,j]`and whose coefficient field `K` is by default `QQ`.

If `cached` is `true`, the internally generated polynomial ring will be cached.

## Examples

```jldoctest
julia> R = gaussian_ring(3)
Gaussian ring over Rational field in 6 variables
s[1, 1], s[1, 2], s[1, 3], s[2, 2], s[2, 3], s[3, 3]
```
"""
function gaussian_ring(n::Int; s_var_name::VarName="s", K::Field=QQ, cached=false)
  varindices = [Tuple([i,j]) for i in 1:n for j in i:n]
  varnames = ["$(s_var_name)[$(i), $(j)]" for (i,j) in varindices]
  S, s = polynomial_ring(K, varnames; cached=cached)
  d = Dict([varindices[i] => s[i] for i in 1:length(varindices)])
  cov_matrix = matrix([[i < j ? d[i,j] : d[j,i] for j in 1:n] for i in 1:n])
  GaussianRing(S, d, cov_matrix)
end

function Base.show(io::IO, R::GaussianRing)
  coeffs = base_ring(R.ring)
  k = ngens(R.ring)
  print(io, "Gaussian ring over $(coeffs) in $(k) variables", "\n", chop(string(gens(R.ring)), head = 16, tail = 1))
end

# functions for getting attributes of Gaussian rings

@doc raw"""
    ring(R::GaussianRing)

Return the multivariate polynomial ring inside `R`.

## Examples

```jldoctest
julia> R = gaussian_ring(3)
Gaussian ring over Rational field in 6 variables
s[1, 1], s[1, 2], s[1, 3], s[2, 2], s[2, 3], s[3, 3]

julia> ring(R)
Multivariate polynomial ring in 6 variables s[1, 1], s[1, 2], s[1, 3], s[2, 2], ..., s[3, 3]
  over rational field
```
"""
function ring(R::GaussianRing)
  R.ring
end

@doc raw"""
    gens(R::GaussianRing)

Return the generators of the multivariate polynomial ring inside the GaussianRing as a dictionary which can be easily indexed

## Examples

```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> R = gaussian_ring(3)
Gaussian ring over Rational field in 6 variables
s[1, 1], s[1, 2], s[1, 3], s[2, 2], s[2, 3], s[3, 3]

julia> gens(R)
Dict{Tuple{Int64, Int64}, QQMPolyRingElem} with 6 entries:
  (1, 2) => s[1, 2]
  (1, 1) => s[1, 1]
  (3, 3) => s[3, 3]
  (1, 3) => s[1, 3]
  (2, 2) => s[2, 2]
  (2, 3) => s[2, 3]
```
"""
function gens(R::GaussianRing)
  R.gens
end

@doc raw"""
    covariance_matrix(R::GaussianRing)

Return the covariance matrix associated to `R` as a matrix over the underlying polynomial ring of `R`

## Examples

```jldoctest
julia> R = gaussian_ring(3)
Gaussian ring over Rational field in 6 variables
s[1, 1], s[1, 2], s[1, 3], s[2, 2], s[2, 3], s[3, 3]

julia> covariance_matrix(R)
[s[1, 1]   s[1, 2]   s[1, 3]]
[s[1, 2]   s[2, 2]   s[2, 3]]
[s[1, 3]   s[2, 3]   s[3, 3]]
```
"""
function covariance_matrix(R::GaussianRing)
  R.covariance_matrix
end

###################################################################################
#
#       Directed Gaussian Graphical Models
#
###################################################################################

@doc raw"""
    graphical_model(G::Graph{Directed}, S::GaussianRing; l_var_name::VarName="l", w_var_name::VarName="w", cached=false)

A parametric statistical model associated to a directed acyclic graph.
It contains a directed acylic graph `G`, a GaussianRing `S` where the vanishing ideal of the model naturally lives, 
and a parameter ring whose variables `l[i,j]` correspond to the directed edges `i -> j` in `G`. 

If `cached` is `true`, the internally generated polynomial ring will be cached.

## Examples

```jldoctest
julia> M = graphical_model(graph_from_edges(Directed, [[1,2], [2,3]]), gaussian_ring(3))
Gaussian graphical model on a directed graph with edges:
(1, 2), (2, 3)
```
"""
function graphical_model(G::Graph{Directed}, S::GaussianRing; l_var_name::VarName="l", w_var_name::VarName="w", cached=false)
  l_indices = [Tuple([src(e),dst(e)]) for e in edges(G)]
  w_indices = vertices(G)
  R, l, w = polynomial_ring(QQ, ["$(l_var_name)[$(i), $(j)]" for (i,j) in l_indices], ["$(w_var_name)[$(v)]" for v in w_indices]; cached=cached)
  d = Dict([l_indices[i] => l[i] for i in 1:length(l_indices)])
  GraphicalModel(G, S, R, [d, w])
end

function Base.show(io::IO, M::GraphicalModel{Graph{Directed}, GaussianRing})
  G = M.graph
  E = [(src(e), dst(e)) for e in edges(G)]

  if length(E) == 0
    print(io, "Gaussian graphical model on the empty directed graph")
  end

  if length(E) > 0
    print(io, "Gaussian graphical model on a directed graph with edges:", "\n", string(E)[2:end-1])
  end
end

@doc raw"""
    directed_edges_matrix(M::GraphicalModel{Graph{Directed}, GaussianRing})

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
function directed_edges_matrix(M::GraphicalModel{Graph{Directed}, GaussianRing})
  G = graph(M)
  l = param_gens(M)[1]
  L = matrix(M.param_ring, [[has_edge(G, i, j) ? l[i,j] : 0 for j in vertices(G)] for i in vertices(G)])
end

@doc raw"""
    error_covariance_matrix(M::GraphicalModel{Graph{Directed}, GaussianRing})

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
function error_covariance_matrix(M::GraphicalModel{Graph{Directed}, GaussianRing})
  G = graph(M)
  w = param_gens(M)[2]
  W = matrix(M.param_ring, [[i == j ? w[i] : 0 for j in vertices(G)] for i in vertices(G)])
end

@doc raw"""
    parametrization(M::GraphicalModel{Graph{Directed}, GaussianRing})

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
function parametrization(M::GraphicalModel{Graph{Directed}, GaussianRing})
  S = ring(M)
  R = param_ring(M)
  G = graph(M)
  L = directed_edges_matrix(M)
  Id = identity_matrix(R, n_vertices(G))
  W = error_covariance_matrix(M)
  Sigma = transpose(inv(Id-L))*W*inv(Id-L)
  hom(S.ring, R, reduce(vcat, [[Sigma[i,j] for j in i:n_vertices(G)] for i in 1:n_vertices(G)]))
end

###################################################################################
#
#       Undirected Gaussian Graphical Models
#
###################################################################################

@doc raw"""
    graphical_model(G::Graph{Undirected}, S::GaussianRing; k_var_name::VarName="k", cached=false)

A parametric statistical model associated to an undirected graph.
It contains an undirected graph `G`, a GaussianRing `S` where the vanishing ideal of the model naturally lives, 
and a parameter ring whose variables `k[i,j]` correspond to the edges `i - j` in `G`. 

If `cached` is `true`, the internally generated polynomial ring will be cached.

## Examples

```jldoctest
julia> M = graphical_model(graph_from_edges([[1,2], [2,3]]), gaussian_ring(3))
Gaussian graphical model on an undirected graph with edges:
(1, 2), (2, 3)
```
"""
function graphical_model(G::Graph{Undirected}, S::GaussianRing; k_var_name::VarName="k", cached=false)
  V = vertices(G)
  varindices = [Tuple([V[i],V[j]]) for i in 1:length(V) for j in i:length(V) if i == j || has_edge(G, V[i], V[j])]
  R, k = polynomial_ring(QQ, ["$(k_var_name)[$(i), $(j)]" for (i,j) in varindices]; cached=cached)
  d = Dict([varindices[i] => k[i] for i in 1:length(varindices)])
  GraphicalModel(G, S, R, d)
end

function Base.show(io::IO, M::GraphicalModel{Graph{Undirected}, GaussianRing})
  G = M.graph
  E = [(dst(e), src(e)) for e in edges(G)]

  if length(E) == 0
    print(io, "Gaussian graphical model on the empty undirected graph")
  end

  if length(E) > 0
    print(io, "Gaussian graphical model on an undirected graph with edges:", "\n", string(E)[2:end-1])
  end
end

@doc raw"""
    concentration_matrix(M::GraphicalModel{Graph{Undirected}, GaussianRing})

Create the concentration matrix `K` of an undirected Gaussian graphical model which is a symmetric positive definite matrix
whose nonzero entries correspond to the edges of the associated graph.

## Examples

```jldoctest
julia> M = graphical_model(graph_from_edges([[1,2], [2,3]]), gaussian_ring(3))
Gaussian graphical model on an undirected graph with edges:
(1, 2), (2, 3)

julia> concentration_matrix(M)
[k[1, 1]   k[1, 2]         0]
[k[1, 2]   k[2, 2]   k[2, 3]]
[      0   k[2, 3]   k[3, 3]]
```
"""
function concentration_matrix(M::GraphicalModel{Graph{Undirected}, GaussianRing})
  G = graph(M)
  k = param_gens(M)
  K = zero_matrix(param_ring(M), n_vertices(G), n_vertices(G))

  for e in edges(G)
    K[dst(e), src(e)] = k[dst(e), src(e)]
    K[src(e), dst(e)] = k[dst(e), src(e)]
  end

  for v in vertices(G)
    K[v, v] = k[v,v]
  end

  K
end

@doc raw"""
    parametrization(M::GraphicalModel{Graph{Undirected}, GaussianRing})

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
function parametrization(M::GraphicalModel{Graph{Undirected}, GaussianRing})
  G = graph(M)
  S = ring(M)
  R = param_ring(M)
  K = concentration_matrix(M)
  adj = adjugate(K)
  hom(ring(S), fraction_field(R), reduce(vcat, [[adj[i,j]//det(K) for j in i:n_vertices(G)] for i in 1:n_vertices(G)]))
end

@doc raw"""
    vanishing_ideal(M::GraphicalModel{Graph{Undirected}, GaussianRing})

Compute the vanishing ideal of the undirected Gaussian graphical model `M`.
This is done by saturating the ideal given by $ \Sigma K - Id $ by the determinant of $ K $
and then eliminating all variables `k[i,j]` where $ K =$ `concentration_matrix(M)`.

## Examples

```jldoctest
julia> M = graphical_model(graph_from_edges([[1,2], [2,3]]), gaussian_ring(3))
Gaussian graphical model on an undirected graph with edges:
(1, 2), (2, 3)

julia> vanishing_ideal(M)
Ideal generated by
  -s[1, 2]*s[2, 3] + s[1, 3]*s[2, 2]
```
"""
function vanishing_ideal(M::GraphicalModel{Graph{Undirected}, GaussianRing})
  G = graph(M)
  S = ring(M)
  R = param_ring(M)

  # there simply must be a better way to do this but I cannot find one currently
  elim_ring, elim_gens = polynomial_ring(coefficient_ring(R), vcat([string(x) for x in gens(R)], [string(x) for x in gens(ring(S))]))
  inject_R = hom(R, elim_ring, elim_gens[1:ngens(R)])
  inject_S = hom(S.ring, elim_ring, elim_gens[ngens(R)+1:ngens(elim_ring)])

  Sigma = inject_S.(covariance_matrix(S))
  K = inject_R.(concentration_matrix(M))

  elim_ideal = ideal(reduce(vcat, Sigma*K - identity_matrix(elim_ring, n_vertices(G))))
  eliminate(saturation(elim_ideal, ideal([det(K)])), elim_gens[1:ngens(R)])
end
