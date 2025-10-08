###################################################################################
#
#       (Un)Directed Gaussian Graphical Models
#
###################################################################################

@attributes mutable struct GaussianGraphicalModel{T, L} <: GraphicalModel{T, L}
  graph::T
  labelings::L
  varnames::Dict{Symbol,VarName}
  function GaussianGraphicalModel(G::Graph{S},
                                  varnames::Dict{Symbol, VarName}) where S <: GraphTypes
    graph_maps = NamedTuple(_graph_maps(G))
    graph_maps = isempty(graph_maps) ? nothing : graph_maps
    return new{Graph{S}, typeof(graph_maps)}(G, graph_maps, varnames)
  end

  function GaussianGraphicalModel(G::MixedGraph, varnames::Dict{Symbol, VarName})
    #TODO figure out how to deal with labelings on MixedGraphs
    # for now just use Nothing
    return new{MixedGraph, Nothing}(G, nothing, varnames)
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
function gaussian_graphical_model(G::Graph{T}; s_varname::VarName="s", l_varname::VarName="l", w_varname::VarName="w") where T <: Directed
  @req is_acyclic(G) "$G must be acyclic"
  GaussianGraphicalModel(G, Dict{Symbol, VarName}(:s => s_varname, :l => l_varname, :w => w_varname))
end

function gaussian_graphical_model(G::Graph{Undirected}; s_varname::VarName="s", k_varname::VarName="k")
  GaussianGraphicalModel(G, Dict{Symbol, VarName}(:s => s_varname, :k => k_varname))
end

function gaussian_graphical_model(G::MixedGraph; s_varname::VarName="s", l_varname::VarName="l", w_varname::VarName="w")
  GaussianGraphicalModel(G, Dict{Symbol, VarName}(:s => s_varname, :l => l_varname, :w => w_varname))
end

varnames(GM::GaussianGraphicalModel) = GM.varnames

function Base.show(io::IO, M::GaussianGraphicalModel{T, L}) where {T, L}
  io = pretty(io)
  if is_terse(io)
    print(io, "Gaussian Graphical Model on a")
    !(L == Nothing) && print(io, " labelled")
    print(io, " $T")
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
  GraphGenDict
} function parameter_ring(GM::GaussianGraphicalModel{Graph{Directed}, T}; cached=false) where T
  G = graph(GM)
  gen_names = (["$(varnames(GM)[:l])[$(src(e)), $(dst(e))]" for e in edges(G)],
               ["$(varnames(GM)[:w])[$(v)]" for v in vertices(G)])
  R, e_gens, v_gens = polynomial_ring(QQ, gen_names; cached=cached)
  gens_dict = merge(Dict(e => e_gens[i] for (i, e) in enumerate(edges(G))),
                    Dict(v => v_gens[v] for v in 1:n_vertices(G)))
  return R, Dict{Union{Int, Edge}, MPolyRingElem}(gens_dict)
end

@attr Tuple{
  MPolyRing,
  GraphGenDict
} function parameter_ring(GM::GaussianGraphicalModel{MixedGraph, T}; cached=false) where T
  G = graph(GM)
  gen_names = (["$(varnames(GM)[:l])[$(src(e)), $(dst(e))]" for e in edges(G, Directed)],
               ["$(varnames(GM)[:w])[$(v)]" for v in vertices(G)], 
               ["$(varnames(GM)[:w])[$(src(e)), $(dst(e))]" for e in edges(G, Undirected)])
  R, d_gens, v_gens, u_gens = polynomial_ring(QQ, gen_names; cached=cached)
  gens_dict = merge(Dict(e => d_gens[i] for (i, e) in enumerate(edges(G, Directed))),
                    Dict(v => v_gens[v] for v in 1:n_vertices(G)),
                    Dict(e => d_gens[i] for (i, e) in enumerate(edges(G, Undirected))))
  return R, Dict{Union{Int, Edge}, MPolyRingElem}(gens_dict)
end

@doc raw"""
    directed_edges_matrix(M::GaussianGraphicalModel{Graph{T}, L}) where {T <: Union{Directed, Mixed}, L}

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
function directed_edges_matrix(M::GaussianGraphicalModel{<:AbstractGraph{T}}) where {T <: Union{Mixed, Directed}}
  G = graph(M)
  R, gens_dict = parameter_ring(M)
  n =  n_vertices(G)
  lambda = zero_matrix(R, n, n)
  for e in edges(G, Directed)
    lambda[src(e), dst(e)] = gens_dict[e]
  end
  return lambda
end

@doc raw"""
    error_covariance_matrix(M::GaussianGraphicalModel{Graph{Directed}, L}) where L

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
function error_covariance_matrix(M::GaussianGraphicalModel{Graph{Directed}, L}) where L
  G = graph(M)
  R, gens_dict = parameter_ring(M)
  W = diagonal_matrix(R, [gens_dict[i] for i in 1:n_vertices(G)])
end


@doc raw"""
    error_covariance_matrix(M::GaussianGraphicalModel{Graphh{Directed}, L}) where L

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
function error_covariance_matrix(M::GaussianGraphicalModel{MixedGraph, L}) where L
  G = graph(M)
  R, gens_dict = parameter_ring(M)
  W = diagonal_matrix(R, [gens_dict[i] for i in 1:n_vertices(G)])

  for e in edges(G, Undirected)
    # TODO maybe using an ordering would be better so we only need one pair?
    W[dst(e), src(e)] = gens_dict[e]
    W[src(e), dst(e)] = gens_dict[e]
  end

  return W
end

function principal_minor(A::Generic.MatSpaceElem, K::Vector{Int})::MPolyRingElem
  det(A[K,K])
end

@doc raw"""
    vanishing_ideal(M::GaussianGraphicalModel{Graph{Directed}, L} where L

The vanishing ideal of a directed acyclic Gaussian graphical model can be computed
by saturating the conditional independence ideal of its ordered pairwise Markov
property at the principal minors corresponding to all parent sets in the graph.
This saturation-based approach usually outperforms the elimination-based default
for general graphical models.

If the graph is cyclic, we fall back to elimination.
"""
function vanishing_ideal(M::GaussianGraphicalModel{Graph{Directed}, L}) where L
  if !is_acyclic(graph(M))
    return vanishing_ideal(M; algorithm=:eliminate)
  end
  G = graph(M)
  S, _ = model_ring(M)
  A = covariance_matrix(M)
  U = foldr(product,
            [powers_of_element(principal_minor(A, parents(G, v))) for v in vertices(G)];
            init=powers_of_element(S(1)))
  # An especially small Markov property
  T = topological_sort(G)
  J = ideal(S, [ci_polynomial(A, ci_stmt(T[i], T[j], parents(G, T[j])))
                for i in 1:length(T) for j in (i+1):length(T)
                if !has_edge(G, T[i], T[j])])
  loc, iota = localization(S, U)
  saturated_ideal(iota(J))
end

@doc raw"""
    parametrization(M::GaussianGraphicalModel{Graph{Directed}, L}) where L

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
function parametrization(M::GaussianGraphicalModel{<: AbstractGraph{T}, L}) where {T  <: Union{Directed, Mixed}, L}
  S, _ = model_ring(M)
  R, _ = parameter_ring(M)
  G = graph(M)
  lambda = directed_edges_matrix(M)
  Id = identity_matrix(R, n_vertices(G))
  W = error_covariance_matrix(M)
  Sigma = transpose(inv(Id - lambda)) * W * inv(Id - lambda)

  gens_map = reduce(vcat, [[Sigma[i,j] for j in i:n_vertices(G)] for i in 1:n_vertices(G)])
  return hom(S, R, gens_map)
end

@doc raw"""
    ci_structure(M::GaussianGraphicalModel{Graph{Directed}, L})

Return the set of elementary CI statements which are satisfied by every distribution in the
graphical model. This is the same as the `global_markov` property of the underlying graph.
"""
function ci_structure(M::GaussianGraphicalModel{Graph{Directed}, L}) where L
  global_markov(graph(M))
end

###################################################################################
#
#       Undirected parametrization
#
###################################################################################
@attr Tuple{
  MPolyRing,
  GraphGenDict
} function parameter_ring(GM::GaussianGraphicalModel{Graph{Undirected}, T}; cached=false) where T
  G = graph(GM)
  gen_names = (["$(varnames(GM)[:k])[$(src(e)), $(dst(e))]" for e in edges(G)],
              ["$(varnames(GM)[:k])[$(v), $(v)]" for v in vertices(G)])
  R, e_gens, v_gens = polynomial_ring(QQ, gen_names; cached=cached)
  gens_dict = merge(Dict(e => e_gens[i] for (i, e) in enumerate(edges(G))),
                    Dict(v => v_gens[v] for v in 1:n_vertices(G)))
  return R, Dict{Union{Int, Edge}, MPolyRingElem}(gens_dict)
end

@doc raw"""
    concentration_matrix(M::GaussianGraphicalModel{Graph{Undirected}, T}) where T

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
function concentration_matrix(M::GaussianGraphicalModel{Graph{Undirected}, T}) where T
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
    vanishing_ideal(M::GaussianGraphicalModel{Graph{Undirected}, L} where L

It is a well-known theorem that the vanishing ideal is the unique prime ideal
above the pairwise conditional independence ideal which does not contain the
determinant of the covariance matrix. This translates into a saturation problem
which can usually be solved faster than the generic elimination-based approach.
"""
function vanishing_ideal(M::GaussianGraphicalModel{Graph{Undirected}, L}) where L
  G = graph(M)
  S, _ = model_ring(M)
  A = covariance_matrix(M)
  U = powers_of_element(det(A))
  J = ideal(S, [ci_polynomial(A, stmt) for stmt in pairwise_markov(G)])
  loc, iota = localization(S, U)
  saturated_ideal(iota(J))
end

@doc raw"""
    parametrization(M::GaussianGraphicalModel{Graph{Undirected}, Nothing})

Create the polynomial map which parametrizes the undirected Gaussian graphical model `M`. Its parameter space is the set of
$K = $ `concentration_matrix(M)` and the map is matrix inversion $K \mapsto K^{-1}$. The constructed ring map writes out the
entries of $K^{-1}$ using the cofactor formula. The target ring is a localization of the polynomial ring at the determinant
of $K$.

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
function parametrization(M::GaussianGraphicalModel{Graph{Undirected}, T}) where T
  G = graph(M)
  S, _ = model_ring(M)
  R, _ = parameter_ring(M)
  K = concentration_matrix(M)
  Rloc, iota = localization(R, powers_of_element(det(K)))
  adj = adjugate(K)
  hom(S, Rloc, reduce(vcat, [[iota(adj[i,j])/iota(det(K)) for j in i:n_vertices(G)] for i in 1:n_vertices(G)]))
end

@doc raw"""
    ci_structure(M::GaussianGraphicalModel{Graph{Undirected}, T})

Return the set of elementary CI statements which are satisfied by every distribution in the
graphical model. This is the same as the `global_markov` property of the underlying graph.
"""
function ci_structure(M::GaussianGraphicalModel{Graph{Undirected}, T}) where T
  global_markov(graph(M))
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

covariance_matrix(R::GaussianRing) = R.covariance_matrix

function Base.show(io::IO, R::GaussianRing)
  coeffs = base_ring(R.ring)
  k = ngens(R.ring)
  print(io, "$(typeof(R)) over $(coeffs) in $(k) variables", "\n", chop(string(gens(R.ring)), head = 16, tail = 1))
end

@doc raw"""
    random_variables(R::GaussianRing)

Return the vector of `[1, ..., n]` where `n` is the number of random variables in the GaussianRing
(equivalently the `n` for which the covariance matrix is an `n`-by-`n` matrix).

## Examples

```jldoctest
julia> R = gaussian_ring(4)
GaussianRing over Rational field in 10 variables
s[1, 1], s[1, 2], s[1, 3], s[1, 4], s[2, 2], s[2, 3], s[2, 4], s[3, 3], s[3, 4], s[4, 4]

julia> random_variables(R)
4-element Vector{Int64}:
 1
 2
 3
 4
```
"""
function random_variables(R::GaussianRing)
  return collect(1:nrows(covariance_matrix(R)))
end

@doc raw"""
    ci_statements(R::GaussianRing)

Return all the `CIStmt` objects which can be formed on the `random_variables(R)`.

```jldoctest
julia> R = gaussian_ring(3)
Gaussian ring over Rational field in 6 variables
s[1, 1], s[1, 2], s[1, 3], s[2, 2], s[2, 3], s[3, 3]

julia> ci_statements(R)
6-element Vector{CIStmt}:
 [1 _||_ 2 | {}]
 [1 _||_ 2 | 3]
 [1 _||_ 3 | {}]
 [1 _||_ 3 | 2]
 [2 _||_ 3 | {}]
 [2 _||_ 3 | 1]
```
"""
ci_statements(R::GaussianRing) = ci_statements(random_variables(R))

@doc raw"""
    ci_ideal(R::GaussianRing, stmts::Vector{CIStmt})

Return the ideal for the conditional independence statements given by `stmts`.

## Examples

```jldoctest
julia> R = gaussian_ring(3)
Gaussian ring over Rational field in 6 variables
s[1, 1], s[1, 2], s[1, 3], s[2, 2], s[2, 3], s[3, 3]

julia> ci_ideal(R, [CI"1,2|", CI"1,2|3"])
Ideal generated by
  s[1, 2]
  s[1, 2]*s[3, 3] - s[1, 3]*s[2, 3]
```
"""
function ci_ideal(R::GaussianRing, stmts::Vector{CIStmt})
  eqs = Any[]
  A = covariance_matrix(R)
  for stmt in stmts
    k = length(stmt.K)
    B = A[vcat(stmt.I, stmt.K), vcat(stmt.J, stmt.K)]
    append!(eqs, minors(B, k+1))
  end
  return ideal(R.ring, eqs)
end

@doc raw"""
    ci_polynomial(R::GaussianRing, stmt::CIStmt)
    ci_polynomial(A::Generic.MatSpaceElem, stmt::CIStmt)

Return the subdeterminant of the matrix `A` which corresponds to the
elementary CI statement `stmt`. If a `GaussianRing` is given, its
`covariance_matrix` is used.

## Examples

```jldoctest
julia> R = gaussian_ring(3)
Gaussian ring over Rational field in 6 variables
s[1, 1], s[1, 2], s[1, 3], s[2, 2], s[2, 3], s[3, 3]

julia> ci_polynomial(R, [CI"1,2|3"])
s[1, 2]*s[3, 3] - s[1, 3]*s[2, 3]
```
"""
function ci_polynomial(A::Generic.MatSpaceElem, stmt::CIStmt)
  @req is_elementary(stmt) "only elementary CI statements give polynomials - for non-elementary statements use ci_ideal"
  det(A[vcat(stmt.I, stmt.K), vcat(stmt.J, stmt.K)])
end

ci_polynomial(R::GaussianRing, stmt::CIStmt) = ci_polynomial(covariance_matrix(R), stmt)

@doc raw"""
    ci_structure(R::GaussianRing; strategy::Symbol)

Return the set of elementary CI statements which are satisfied by every distribution in the
given Gaussian model. For this method to work, the ring needs to have one of the following
methods:

1. `parametrization(R)` returning a parametrized covariance matrix. The CI structure is determined by evaluating its subdeterminants.
2. `vanishing_ideal(R)` returning the model's vanishing ideal. In this case, the subdeterminant for each CI statement is checked for contained in the vanishing ideal.

The optional argument `strategy` can be used to give preference to one of the above methods using the values `:parametrization` or `:ideal`, respectively.
"""
function ci_structure(R::GaussianRing; strategy::Symbol)
  # TODO: Get rid of hasmethod probing.
  if strategy == :parametrization || hasmethod(parametrization, Tuple{typeof(R)})
    A = parametrization(R)
    return [stmt for stmt in ci_statements(R) if ci_polynomial(A, stmt) == 0]
  elseif strategy == :ideal || hasmethod(vanishing_ideal, Tuple{typeof(R)})
    A = covariance_matrix(R)
    I = vanishing_ideal(R)
    return [stmt for stmt in ci_statements(R) if ci_polynomial(A, stmt) in I]
  end
  error("no method available to compute the CI structure")
end
