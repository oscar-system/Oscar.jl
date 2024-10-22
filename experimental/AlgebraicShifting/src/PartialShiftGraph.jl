"""
  Discovers the nodes of partial shift graph starting from `K`.
"""
function partial_shift_graph_nodes(F::Field, K::SimplicialComplex, W::Vector{WeylGroupElem};)
  current_node = K
  visited_nodes = Set([current_node])
  
  function adjacent_shifts(K::SimplicialComplex)
    shifts_of_K = [exterior_shift(F, K, w) for w in W]
    adjacent_nodes = Set(filter(x -> !isequal(K, x), shifts_of_K))
    return shifts_of_K
  end
  
  unvisited_nodes = adjacent_shifts(current_node)

  while(!isempty(unvisited_nodes))
    current_node = pop!(unvisited_nodes)
    push!(visited_nodes, current_node)
    shifts = adjacent_shifts(current_node)
    # dont visit things twice
    unvisited_nodes = union(unvisited_nodes, setdiff(shifts, visited_nodes))
  end

  return sort(collect(visited_nodes))
end


""" Compute the multi edges, that is, for each complex K  compute which
    other complexes can be reached by applying the partial shifts `deltas`
    to K, and store the shift matrices that give rise to each edge."""
function multi_edges(F::Field,
                     W::WeylGroup,
                     complexes::Vector{Tuple{Int,T}},
                     complex_labels::Dict{T, Int}
) :: Dict{Tuple{Int, Int}, Vector{WeylGroupElem}} where T <: ComplexOrHypergraph
  # For each complex K with index i, compute the shifted complex Δ of K by w for each w ∈ W.
  # For nontrivial Δ, place (i, Δ) → w in a singleton dictionary, and eventually merge all dictionaries
  # to obtain a dictionary (i, Δ) → [w that yield the shift K → Δ]
  reduce(
   (d1, d2) -> mergewith!(vcat, d1, d2), 
    (
      Dict((i, complex_labels[Δ]) => [w])
      for (i, K) in complexes
      for (w, Δ) in ((w, exterior_shift2(K, w; field=F)) for w in W)
      if Δ != K
    );
    init=Dict{Tuple{Int, Int}, Vector{WeylGroupElem}}()
  )
end


"""    construct_full_graph(complexes :: Vector{T} where T <: ComplexOrHypergraph; ...) :: Tuple{Graph{Directed}, EdgeLabels}

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
function partial_shift_graph(F::Field, complexes :: Vector{T} where T <: ComplexOrHypergraph;
                             parallel :: Bool = false, W = nothing) :: Tuple{Graph{Directed}, EdgeLabels}
  # Deal with trivial case
  if length(complexes) <= 1
    return (graph_from_adjacency_matrix(Directed, zeros(length(complexes),length(complexes))), EdgeLabels())
  end

  ns_vertices = Set(n_vertices.(complexes))
  @req length(ns_vertices) == 1 "All complexes are required to have the same number of vertices."
  n = collect(ns_vertices)[1]

  # inverse lookup K → index of K in complexes
  complex_labels = Dict(K => index for (index, K) in enumerate(complexes))

  # set the weyl group to be used to construct the shifts
  if isnothing(W)
     W = weyl_group(:A, n - 1);
  else
    W2 = only(unique(parent.(W)))
    T = root_system_type(root_system(W2))
    @req T[1][1] == :A && T[1][2] == n - 1 "The Weyl group should be of type A_$(n-1) and not $(T[1])."
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
  edge_labels = reduce(
    (d1, d2) -> mergewith!(vcat, d1, d2),
    # a progress bar here would be nice for the users, but this requires adding a dependency to Oscar
    # @showprogress map_function(
    map_function(
      Ks -> multi_edges(F, W, Ks, complex_labels), 
      Iterators.partition(enumerate(complexes), task_size)
    )
  )  

  graph = graph_from_edges(Directed, [[i,j] for (i,j) in keys(edge_labels)])
  return (graph, edge_labels)
end

