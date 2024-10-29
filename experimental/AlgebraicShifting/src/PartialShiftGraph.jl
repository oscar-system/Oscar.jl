const EdgeLabels = Dict{Tuple{Int, Int}, Vector{WeylGroupElem}}

function isless_lex(S1::Set{Set{Int}}, S2::Set{Set{Int}})
  S_diff = collect(symdiff(S1, S2))
  isempty(S_diff) && return false
  set_cmp(a, b) = min(symdiff(a, b)...) in a

  return sort(S_diff;lt=set_cmp)[1] in S1
end

"""
  Given two simplicial complexes `K1`, `K2` return true if
  `K1` is lexicographically less than `K2`
"""
isless_lex(K1::SimplicialComplex, K2::SimplicialComplex) = isless_lex(Set(facets(K1)), Set(facets(K2)))

"""
  Discovers the nodes of partial shift graph starting from `K`.
"""
function partial_shift_graph_vertices(F::Field,
                                      K::SimplicialComplex,
                                      W::Union{WeylGroup, Vector{WeylGroupElem}};)
  current = K
  visited = [current]
  # by properties of algebraic shifting
  # we know that K we be the last in this list since it is sorted
  unvisited = unique(
    x -> Set(facets(x)),
    sort([exterior_shift(F, K, w) for w in W]; lt=isless_lex))[1:end - 1]

  while(!isempty(unvisited))
    current = pop!(unvisited)
    push!(visited, current)
    shifts = unique(
      x -> Set(facets(x)),
      sort([exterior_shift(F, current, w) for w in W]; lt=isless_lex))[1:end - 1]

    # dont visit things twice
    new_facets = filter(x -> !(x in Set.(facets.(visited))), Set.(facets.(shifts)))
    unvisited = simplicial_complex.(
      union(Set.(facets.(unvisited)), new_facets))
  end
  return sort(collect(visited); lt=isless_lex)
end

""" Compute the multi edges, that is, for each complex `K`  compute which
    other complexes can be reached by applying the partial shifts `deltas`
    to K, and store the shift matrices that give rise to each edge."""
function multi_edges(F::Field,
                     W::WeylGroup,
                     complexes::Vector{Tuple{Int,T}},
                     complex_labels::Dict{Set{Set{Int}}, Int}
) :: Dict{Tuple{Int, Int}, Vector{WeylGroupElem}}  where T <: ComplexOrHypergraph;
  # For each complex K with index i, compute the shifted complex delta of K by w for each w ∈ W.
  # For nontrivial delta, place (i, delta) → w in a singleton dictionary, and eventually merge all dictionaries
  # to obtain a dictionary (i, delta) → [w that yield the shift K → delta]
  reduce(
  (d1, d2) -> mergewith!(vcat, d1, d2), 
  (
    Dict((i, complex_labels[Set(facets(delta))]) => [w])
    for (i, K) in complexes
      for (w, delta) in ((w, exterior_shift(F, K, w)) for w in W)
        if Set(facets(delta)) != Set(facets(K)));
    init=Dict{Tuple{Int, Int}, Vector{WeylGroupElem}}()
  )
end

@doc raw"""
     partial_shift_graph(complexes::Vector{Simplicialcomplex})
     partial_shift_graph(complexes::Vector{Uniformhypergraph})

Constructs the partial shift graph on `complexes`.

Returns a tuple `(G, D)`, where `G` is a directed graph whose vertices correspond to the `complexes`,
such that there is an edge `K → L` if there exists a shift matrix `w` such that `L` is the partial generic shift of `K` by `w`.
If `K` and `K` are the `i`th and `j`th entry of `complexes`, resp., 
`D[i,j]` contains all `w ∈ W` such that `L` is the partial generic shift of `K` by `w`.

# Arguments
- `complexes`: A vector of simplicial complexes of uniform hypergraphs (all belonging to the same ``Γ(n, k, l)``).
- `parallel :: Bool` (default: `false`) run the process in parrallel using the `Distributed` package; the setup is left to the user.
- `W`: The user may provide a list `W` of Weyl group elements to be used to construct the shifts.
  `W` must be a subset the (same instance of the) symmetric group of the same order as the complexes.
  If `W` is not provided, the function will use the symmetric group of the same order as the complexes.

# Examples
```
julia> Γ(n,k,l) = uniform_hypergraph.(subsets(subsets(n, k), l), n)
julia> Ks = Γ(4,2,5)
julia> G, D = construct_full_graph(Ks)
```
"""
function partial_shift_graph(F::Field, complexes::Vector{T}, W::Union{Nothing, WeylGroup, Vector{WeylGroupElem}} = nothing;
                             parallel::Bool = false) :: Tuple{Graph{Directed}, EdgeLabels}  where T <: ComplexOrHypergraph;
  # Deal with trivial case
  if length(complexes) <= 1
    return (graph_from_adjacency_matrix(Directed, zeros(length(complexes),length(complexes))), EdgeLabels())
  end

  ns_vertices = Set(n_vertices.(complexes))
  @req length(ns_vertices) == 1 "All complexes are required to have the same number of vertices."
  n = collect(ns_vertices)[1]

  # inverse lookup K → index of K in complexes
  complex_labels = Dict(Set(facets(K)) => index for (index, K) in enumerate(complexes))

  # set the weyl group to be used to construct the shifts
  if isnothing(W)
     W = weyl_group(:A, n - 1);
  else
    W2 = only(unique(parent.(W)))
    rs_type = root_system_type(root_system(W2))
    @req rs_type[1][1] == :A && rs_type[1][2] == n - 1 "Only Weyl groups type A_$(n-1) are currently support and received type $(T[1])."
  end
  
  task_size = 1
  map_function = map
  if parallel
    # setup parallel parameters
    channels = Oscar.params_channels(Union{RootSystem, WeylGroup, Vector{SimplicialComplex}, Vector{UniformHypergraph}})
    # setup parents needed to be sent to each process
    Oscar.put_params(channels, root_system(W))
    Oscar.put_params(channels, W)
    map_function = pmap
  end
  
  # Basically, the following does the same as the following, but with a progress indicator:
  # edge_tuples = multi_edges(W, collect(enumerate(complexes)), complex_labels)
  edge_labels =  multi_edges(F, W, collect(enumerate(complexes)), complex_labels)
  #edge_labels = reduce(
  #  (d1, d2) -> mergewith!(vcat, d1, d2),
  #  # a progress bar here would be nice for the users, but this requires adding a dependency to Oscar
  #  # @showprogress map_function(
  #  map_function(
  #    Ks -> multi_edges(F, W, Ks, complex_labels), 
  #    Iterators.partition(enumerate(complexes), task_size)
  #  )
  #)  

  graph = graph_from_edges(Directed, [[i,j] for (i,j) in keys(edge_labels)])
  return (graph, edge_labels)
end

