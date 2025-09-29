# TODO: change Vector -> Set
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
isless_lex(K1::ComplexOrHypergraph, K2::ComplexOrHypergraph) = isless_lex(Set(facets(K1)), Set(facets(K2)))

@doc raw"""
    partial_shift_graph_vertices(F::Field,K::SimplicialComplex, W::Union{WeylGroup, Vector{WeylGroupElem}};)
    partial_shift_graph_vertices(F::Field,K::UniformHypergraph, W::Union{WeylGroup, Vector{WeylGroupElem}};)


Given a field `F` find the vertices of the partial shift graph starting from `K`
and discoverable from elements in `W`.
Return a `Vector{SimplicialCompplex}` ordered lexicographically.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> K = simplicial_complex([[1, 2], [2, 3], [3, 4]])
Abstract simplicial complex of dimension 1 on 4 vertices

julia> shifts = partial_shift_graph_vertices(QQ, K, weyl_group(:A, 3))
6-element Vector{SimplicialComplex}:
 Abstract simplicial complex of dimension 1 on 4 vertices
 Abstract simplicial complex of dimension 1 on 4 vertices
 Abstract simplicial complex of dimension 1 on 4 vertices
 Abstract simplicial complex of dimension 1 on 4 vertices
 Abstract simplicial complex of dimension 1 on 4 vertices
 Abstract simplicial complex of dimension 1 on 4 vertices

julia> facets.(shifts)
6-element Vector{Vector{Set{Int64}}}:
 [Set([3, 1]), Set([4, 1]), Set([2, 1])]
 [Set([3, 1]), Set([2, 3]), Set([2, 1]), Set([4])]
 [Set([3, 1]), Set([4, 2]), Set([2, 1])]
 [Set([4, 3]), Set([3, 1]), Set([2, 1])]
 [Set([2, 1]), Set([2, 3]), Set([4, 2])]
 [Set([2, 1]), Set([2, 3]), Set([4, 3])]
```
"""
function partial_shift_graph_vertices(F::Field,
                                      K::SimplicialComplex,
                                      W::Union{WeylGroup, Vector{WeylGroupElem}})
  current = K
  visited = [current]
  phi = isomorphism(PermGroup, parent(first(W)))
  # by properties of algebraic shifting
  # we know that K will be the last in this sorted list
  # sorting here should also speed up unique according to julia docs
  unvisited = unique(
    x -> Set(facets(x)),
    sort([exterior_shift(F, K, phi(w)) for w in W]; lt=isless_lex))[1:end - 1]

  while !isempty(unvisited)
    current = pop!(unvisited)
    push!(visited, current)
    shifts = unique(
      x -> Set(facets(x)),
      sort([exterior_shift(F, current, phi(w)) for w in W]; lt=isless_lex))[1:end - 1]

    # dont visit things twice
    new_facets = filter(x -> !(x in Set.(facets.(visited))), Set.(facets.(shifts)))
    unvisited = simplicial_complex.(
      union(Set.(facets.(unvisited)), new_facets))
  end
  return sort(collect(visited); lt=isless_lex)
end

function partial_shift_graph_vertices(F::Field,
                                      K::UniformHypergraph,
                                      W::Union{WeylGroup, Vector{WeylGroupElem}})
  return partial_shift_graph_vertices(F, simplicial_complex(K), W)
end

""" Compute the multi edges, that is, for each complex `K`  compute which
    other complexes can be reached by applying the partial shifts `deltas`
    to K, and store the shift matrices that give rise to each edge."""
function multi_edges(F::Field,
                     permutations::Vector{PermGroupElem},
                     complexes::Vector{Tuple{Int,T}},
                     complex_labels::Dict{Set{Set{Int}}, Int}
) where T <: ComplexOrHypergraph;
  # For each complex K with index i, compute the shifted complex delta of K by w for each w in W.
  # For nontrivial delta, place (i, delta) → w in a singleton dictionary, and eventually merge all dictionaries
  # to obtain a dictionary (i, delta) → [w that yield the shift K → delta]
  return reduce(
  (d1, d2) -> mergewith!(vcat, d1, d2), (
    Dict((i, complex_labels[Set(facets(delta))]) => [p])
    for (i, K) in complexes
      for (p, delta) in ((p, exterior_shift(F, K, p)) for p in permutations)
        if !issetequal(facets(delta), facets(K)));
    init=Dict{Tuple{Int, Int}, Vector{PermGroupElem}}()
  ) :: Dict{Tuple{Int, Int}, Vector{PermGroupElem}}
end

@doc raw"""
    partial_shift_graph(F::Field, complexes::Vector{Simplicialcomplex}; parallel=false, show_progress=true)
    partial_shift_graph(F::Field, complexes::Vector{Uniformhypergraph}; parallel=false, show_progress=true)
    partial_shift_graph(F::Field, complexes::Vector{Simplicialcomplex}, W::Union{WeylGroup, Vector{WeylGroupElem}}; parallel=false, show_progress=true)
    partial_shift_graph(F::Field, complexes::Vector{Uniformhypergraph}, W::Union{WeylGroup, Vector{WeylGroupElem}}; parallel=false, show_progress=true)

Construct the partial shift graph on `complexes`.

Return a tuple `(G, EL, VL)`, where `G` is a `Graph{Directed}`, `EL` is a `Dict{Tuple{Int Int}, Vector{Weylgroupelem}` and
`VL` is a lexicographically sorted `complexes`, hence is either a `Vector{SimplicialComplex}` or `Vector{Uniformhypergraph}`.
`EL` are the edges labels and `VL` are the vertex labels.
There is an edge from the vertex labelled `K` to the vertex labelled `L` if `L` is the partial shift of `K` by some `w` in `W`.
If `K` and `L` are the `i`th and `j`th entry of `VL`, resp.,
`EL[i,j]` contains all `w` in `W` such that `L` is the partial generic shift of `K` by `w`.

# Arguments
- `complexes`: A vector of simplicial complexes or uniform hypergraphs (all should have the same number of vertices).
- `W`: The user may provide a list `W` of Weyl group elements to be used to construct the shifts.
  All elements of `W` should have the same parent.
  `W` must be a subset the (same instance of the) symmetric group of order equal to the number of vertices of a complex in complexes (they should all be equal).
  If `W` is not provided, the function will use the symmetric group of the same order as vertices in each complex complexes.
- `parallel`: run the process in parrallel using the `Distributed` package; make sure to do `@everywhere using Oscar`.
- `show_progress`: Set to `true` by default, can be used to suppress progress meter.
- `task_size`: While running in parallel serialization might become a bottleneck,
setting a higher task_size should increase performance if the processes are not running at maximum capacity
See [D-VJL24](@cite) for background.


# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> gamma(n,k,l) = uniform_hypergraph.(subsets(subsets(n, k), l), n)
gamma (generic function with 1 method)

julia> Ks = gamma(4,2,5)
6-element Vector{UniformHypergraph}:
 UniformHypergraph(4, 2, [[1, 2], [1, 3], [1, 4], [2, 3], [3, 4]])
 UniformHypergraph(4, 2, [[1, 2], [1, 3], [2, 3], [2, 4], [3, 4]])
 UniformHypergraph(4, 2, [[1, 2], [1, 3], [1, 4], [2, 3], [2, 4]])
 UniformHypergraph(4, 2, [[1, 2], [1, 4], [2, 3], [2, 4], [3, 4]])
 UniformHypergraph(4, 2, [[1, 2], [1, 3], [1, 4], [2, 4], [3, 4]])
 UniformHypergraph(4, 2, [[1, 3], [1, 4], [2, 3], [2, 4], [3, 4]])

julia> G, EL, VL = partial_shift_graph(QQ, Ks; show_progress=false);

julia> collect(edges(G))
14-element Vector{Edge}:
 Edge(2, 1)
 Edge(3, 1)
 Edge(3, 2)
 Edge(4, 1)
 Edge(4, 2)
 Edge(5, 1)
 Edge(5, 2)
 Edge(5, 3)
 Edge(5, 4)
 Edge(6, 1)
 Edge(6, 2)
 Edge(6, 3)
 Edge(6, 4)
 Edge(6, 5)

julia> EL[6, 5]
4-element Vector{WeylGroupElem}:
 s2
 s1 * s2
 s3 * s2
 s1 * s3 * s2

julia> facets.(VL[[6, 5]])
2-element Vector{Vector{Set{Int64}}}:
 [Set([3, 1]), Set([4, 1]), Set([2, 3]), Set([4, 2]), Set([4, 3])]
 [Set([2, 1]), Set([4, 1]), Set([2, 3]), Set([4, 2]), Set([4, 3])]
```
"""
function partial_shift_graph(F::Field, complexes::Vector{T},
                             W::Union{WeylGroup, Vector{WeylGroupElem}};
                             parallel::Bool = false,
                             show_progress::Bool = true,
                             task_size::Int=100) where T <: ComplexOrHypergraph
  # see TODO above about changing EdgeLabels type
  # Deal with trivial case
  if length(complexes) == 1
    @req is_shifted(complexes[1]) "The list of complexes should be closed under shifting by elements of W"
    return (
      graph_from_adjacency_matrix(Directed, zeros(length(complexes),length(complexes))),
      EdgeLabels(),
      complexes) :: Tuple{Graph{Directed}, EdgeLabels, Vector{T}}
  end
  
  # maybe we provide a flag to skip if the complexes are already sorted?
  complexes = sort(complexes;lt=Oscar.isless_lex)

  ns_vertices = Set(n_vertices.(complexes))
  @req length(ns_vertices) == 1 "All complexes are required to have the same number of vertices."
  n = collect(ns_vertices)[1]

  # inverse lookup K → index of K in complexes
  complex_labels = Dict(Set(facets(K)) => index for (index, K) in enumerate(complexes))

  W2 = only(unique(parent.(W)))
  rs_type = root_system_type(root_system(W2))
  @req rs_type[1][1] == :A && rs_type[1][2] == n - 1 "Only Weyl groups type A_$(n-1) are currently support and received type $(T[1])."

  phi = isomorphism(PermGroup, parent(first(W)))
  # oscar runs tests in parallel, which can make pmap run in parallel even if parallel=false
  # next line fixes this issue
  map_function = map
  if parallel
    # setup parallel parameters
    # this should be updated to use the parallel framework at some point?
    channels = [RemoteChannel(()->Channel{Any}(32), i) for i in workers()]
    # setup parents needed to be sent to each process
    map(channel -> put_type_params(channel, codomain(phi)), channels)
    map_function = pmap
  end
  try
    # here to so that this isn't applied on the worker
    # since we cannot serialize Gap maps 
    PG = phi.(W)
    if show_progress
      edge_labels = reduce((d1, d2) -> mergewith!(vcat, d1, d2),
                           @showprogress map_function(
                             Ks -> multi_edges(F, PG, Ks, complex_labels),
                             Iterators.partition(enumerate(complexes), task_size)))
    else
      edge_labels = reduce((d1, d2) -> mergewith!(vcat, d1, d2),
                           map_function(
                             Ks -> multi_edges(F, PG, Ks, complex_labels),
                             Iterators.partition(enumerate(complexes), task_size)))
    end
    graph = graph_from_edges(Directed, [[i,j] for (i,j) in keys(edge_labels)])
    return (graph,
            Dict(k => inv(phi).(v) for (k, v) in edge_labels),
            complexes) :: Tuple{Graph{Directed}, EdgeLabels, Vector{T}}
  catch e
    if e isa KeyError
      error("The list of complexes should be closed under shifting by elements of W")
    end
    rethrow(e)
  end
end

function partial_shift_graph(F::Field, complexes::Vector{T};
                             kwargs...) where T <: ComplexOrHypergraph
  # Deal with trivial case
  if length(complexes) <= 1
    return (graph_from_adjacency_matrix(Directed, zeros(length(complexes),length(complexes))), EdgeLabels())
  end

  n = n_vertices(complexes[1])
  W = weyl_group(:A, n - 1)
  return partial_shift_graph(F, complexes, W; kwargs...)
end

@doc raw"""
    contracted_partial_shift_graph(G::Graph{Directed}, edge_labels::Dict{Tuple{Int, Int}, Vector{WeylGroupElem}})

Return a triple `(CG, S, P)`, where `CG` is a graph that contains a vertex `v` for every vertex `S[v]` in `G`.
`S` is a list of indices for the sinks in the original graph `G`.
A vertex `i` is in `P[s]` if there exists an edge from `i` to `s` in `G` with `w0` in its edge label,
in this way `P` is a partition of the vertices of the orignal graph `G`.
There is an edge from `s` to `t`  in `CG` whenever there is an edge from `i` to `j` in `G` and `i` in `P[s]` and `j` in `P[t]`.
See [D-VJL24](@cite) for background.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> gamma(n,k,l) = uniform_hypergraph.(subsets(subsets(n, k), l), n)
gamma (generic function with 1 method)

julia> Ks = gamma(4,2,5)
6-element Vector{UniformHypergraph}:
 UniformHypergraph(4, 2, [[1, 2], [1, 3], [1, 4], [2, 3], [3, 4]])
 UniformHypergraph(4, 2, [[1, 2], [1, 3], [2, 3], [2, 4], [3, 4]])
 UniformHypergraph(4, 2, [[1, 2], [1, 3], [1, 4], [2, 3], [2, 4]])
 UniformHypergraph(4, 2, [[1, 2], [1, 4], [2, 3], [2, 4], [3, 4]])
 UniformHypergraph(4, 2, [[1, 2], [1, 3], [1, 4], [2, 4], [3, 4]])
 UniformHypergraph(4, 2, [[1, 3], [1, 4], [2, 3], [2, 4], [3, 4]])

julia> G, EL, VL = partial_shift_graph(QQ, Ks; show_progress=false);

julia> contracted_partial_shift_graph(G, EL)
(Directed graph with 1 nodes and 0 edges, [1], [[5, 4, 6, 2, 3, 1]])
```
"""
function contracted_partial_shift_graph(G::Graph{Directed}, edge_labels::Dict{Tuple{Int, Int}, Vector{WeylGroupElem}})
  n = n_vertices(G)
  W = parent(first(first(values(edge_labels))))
  w0 = longest_element(W)

  # all arrows corresponding to edges that contain w0 in their edge labels
  w0_action = Dict(i => j for ((i,j), ws) in edge_labels if w0 in ws)

  sinks = findall(iszero, degree(G))
  sinks_indices = Dict(s => i for (i,s) in enumerate(sinks))
  for i in 1:n
    if !haskey(w0_action, i)
      @req i in sinks "Vertex $i is not a sink, but has no outbound edge with w0 in its edge label."
      w0_action[i] = i
    end
  end

  p = [Int[] for _ in sinks]
  for (i, s) in w0_action
    push!(p[sinks_indices[s]], i)
  end

  return (
    graph_from_edges(Directed, [
      [sinks_indices[s],sinks_indices[t]]
      for (s,t) in (
        (w0_action[i], (haskey(w0_action, j) ? w0_action[j] : j))
        for (i, j) in keys(edge_labels)
      )
      if s != t
    ], length(sinks)),
    sinks,
    p
  )
end
