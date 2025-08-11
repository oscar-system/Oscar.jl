function pm_object(P::PartiallyOrderedSet)
  return P.pm_poset
end

function partially_ordered_set(pp::Polymake.BigObject)
  @req startswith(Polymake.type_name(pp), "PartiallyOrderedSet") "Not a polymake PartiallyOrderedSet: $(Polymake.type_name(pp))"
  return PartiallyOrderedSet(pp)
end

@doc raw"""
    partially_ordered_set(covrels::Union{SMat, AbstractMatrix{<:IntegerUnion}})

Construct a partially ordered set from covering relations `covrels`, given in
the form of the adjacency matrix of a Hasse diagram. The covering relations
must be given in topological order, i.e., `covrels` must be strictly upper
triangular.

The construction only distinguishes zero coefficients (to indicate the absence) and nonzero coefficients (to indicate the presence of a covering relation).

# Examples
```jldoctest
julia> a2_adj_cov = [
           0 1 1 0 0 0
           0 0 0 1 1 0
           0 0 0 1 1 0
           0 0 0 0 0 1
           0 0 0 0 0 1
           0 0 0 0 0 0
         ];

julia> pa2 = partially_ordered_set(a2_adj_cov)
Partially ordered set of rank 3 on 6 elements
```
"""
partially_ordered_set(covrels::AbstractMatrix{<:IntegerUnion}) =
  partially_ordered_set((!iszero).(covrels))

partially_ordered_set(covrels::IncidenceMatrix) =
  partially_ordered_set_from_adjacency(covrels)

partially_ordered_set(covrels::Union{<:AbstractMatrix{Bool},BitMatrix}) =
  partially_ordered_set_from_adjacency(incidence_matrix(covrels))

function partially_ordered_set(sm::SMat)
  incmat = incidence_matrix(nrows(sm), ncols(sm))
  for i in 1:nrows(sm)
    for (k, _) in sm[i]
      incmat[i, k] = 1
    end
  end
  return partially_ordered_set_from_adjacency(incmat)
end

function partially_ordered_set_from_adjacency(covrels::IncidenceMatrix)
  return partially_ordered_set(Graph{Directed}(Polymake.Graph{Directed}(covrels)))
end

@doc raw"""
    partially_ordered_set_from_inclusions(I::IncidenceMatrix)

Construct an inclusion based partially ordered with rows of `I` as co-atoms.
That is, the elements of the poset are the given sets, their intersections, and the union of all rows as greatest element.

# Examples
```jldoctest
julia> im = vertex_indices(facets(simplex(3)))
4×4 IncidenceMatrix
 [1, 3, 4]
 [1, 2, 4]
 [1, 2, 3]
 [2, 3, 4]

julia> partially_ordered_set_from_inclusions(im)
Partially ordered set of rank 4 on 16 elements
```
"""
function partially_ordered_set_from_inclusions(I::IncidenceMatrix)
  # we want to make sure it is always built from the bottom up
  # and just need to specify some arbitrary upper bound for the rank
  ppos = Polymake.call_function(:polytope, :lower_hasse_diagram, I, ncols(I))
  pos = partially_ordered_set(ppos)
  pos.atomlabels = collect(1:ncols(I))
  pos.artificial_top = true
  return pos
end

_pmdec(elem::Int, r::Int) = Polymake.BasicDecoration(Polymake.Set{Int}(elem), r)
_pmdec(r::Int) = Polymake.BasicDecoration(Polymake.Set{Int}(), r)

@doc raw"""
    partially_ordered_set(g::Graph{Directed})

Construct a partially ordered set from a directed graph describing the Hasse diagram.
The graph must be acyclic, edges must be directed from minimal to maximal elements.

If there is no unique minimal (maximal) element in the graph, an artificial least (greatest) element is added to the internal datastructure.

# Examples
```jldoctest
julia> g = graph_from_edges(Directed, [[2,1],[2,4],[1,3],[4,3]])
Directed graph with 4 nodes and the following edges:
(1, 3)(2, 1)(2, 4)(4, 3)

julia> pos = partially_ordered_set(g)
Partially ordered set of rank 2 on 4 elements
```
"""
function partially_ordered_set(gin::Graph{Directed})
  @req n_vertices(gin) >= 2 "Graph must have at least two nodes"
  g = deepcopy(gin)
  dec = Polymake.NodeMap{Directed,Polymake.BasicDecoration}(pm_object(g))
  minimal = findall(iszero, indegree(g))
  maximal = findall(iszero, outdegree(g))
  @req length(minimal) > 0 "no minimal element found"
  @req length(maximal) > 0 "no maximal element found"
  pos = 1
  if length(minimal) > 1
    # need artificial bottom node
    add_vertex!(g)
    least_element = n_vertices(g)
    add_edge!.(Ref(g), Ref(least_element), minimal)
    abottom = true
    maxrank = -1
  else
    least_element = only(minimal)
    abottom = false
    maxrank = 0
  end
  queue = [least_element]
  dec[least_element] = _pmdec(least_element, maxrank)
  # we need to generate a valid rank function
  visited = BitSet(queue)
  # we keep the elements in the queue for permuting the graph at the end
  while pos <= length(queue)
    src = queue[pos]
    r = Polymake.decoration_rank(dec[src])
    for t in outneighbors(g, src)
      rr = max(r + 1, Polymake.decoration_rank(dec[t]))
      maxrank = max(maxrank, rr)
      dec[t] = _pmdec(t, rr)
      if !(t in visited)
        push!(queue, t)
        push!(visited, t)
      end
    end
    pos += 1
  end
  if length(maximal) > 1
    # need artificial top
    add_vertex!(g)
    greatest_element = n_vertices(g)
    add_edge!.(Ref(g), maximal, Ref(greatest_element))
    dec[greatest_element] = _pmdec(greatest_element, maxrank + 1)
    atop = true
    push!(queue, greatest_element)
  else
    atop = false
  end

  Polymake.call_function(
    :common, :permute_graph, pm_object(g), Polymake.to_zero_based_indexing(queue)
  )
  pos = Polymake.graph.PartiallyOrderedSet{Polymake.BasicDecoration}(;
    ADJACENCY=pm_object(g),
    DECORATION=dec,
    TOP_NODE=length(queue) - 1,
    BOTTOM_NODE=0,
  )
  opos = PartiallyOrderedSet(pos)
  opos.artificial_bottom = abottom
  opos.artificial_top = atop
  return opos
end

@doc raw"""
    partially_ordered_set(g::Graph{Directed}, node_ranks::Dict{Int,Int})

Construct a partially ordered set from a directed graph describing the Hasse diagram.
The graph must be acyclic.
The dictionary `node_ranks` must give a valid rank for each node in the graph, strictly
increasing from the minimal elements to maximal elements.
The rank difference between two adjacent nodes may be larger than one.

If there is no unique minimal (maximal) element in the graph, an artificial least (greatest) element is added to the internal datastructure.

# Examples
```jldoctest
julia> g = graph_from_edges(Directed, [[2,1],[2,4],[1,3],[4,3],[3,6],[2,5],[5,6]])
Directed graph with 6 nodes and the following edges:
(1, 3)(2, 1)(2, 4)(2, 5)(3, 6)(4, 3)(5, 6)

julia> pos = partially_ordered_set(g, Dict(2=>0, 1=>1, 4=>2, 3=>3, 5=>2, 6=>4))
Partially ordered set of rank 4 on 6 elements
```
"""
function partially_ordered_set(gin::Graph{Directed}, node_ranks::Dict{Int,Int})
  @req n_vertices(gin) >= 2 "Graph must have at least two nodes"
  g = deepcopy(gin)
  dec = Polymake.NodeMap{Directed,Polymake.BasicDecoration}(pm_object(g))
  for n in 1:n_vertices(g)
    Polymake._set_entry(dec, n - 1, _pmdec(n, node_ranks[n]))
  end
  sorted = sort(collect(node_ranks); by=last)
  sortednodes = first.(sorted)
  if last(sorted[1]) == last(sorted[2])
    # need artificial bottom node
    minimal = first.(filter(x -> last(x) == last(sorted[1]), sorted))
    add_vertex!(g)
    least_element = n_vertices(g)
    add_edge!.(Ref(g), Ref(least_element), minimal)
    dec[least_element] = _pmdec(least_element, -1)
    abottom = true
    pushfirst!(sortednodes, least_element)
  else
    abottom = false
  end
  if last(sorted[end]) == last(sorted[end - 1])
    # need artificial top node
    maximal = first.(filter(x -> last(x) == last(sorted[end]), sorted))
    add_vertex!(g)
    greatest_element = n_vertices(g)
    add_edge!.(Ref(g), maximal, Ref(greatest_element))
    dec[greatest_element] = _pmdec(greatest_element, last(sorted[end]) + 1)
    atop = true
    push!(sortednodes, greatest_element)
  else
    atop = false
  end
  Polymake.call_function(
    :common, :permute_graph, pm_object(g), Polymake.to_zero_based_indexing(sortednodes)
  )
  pos = Polymake.graph.PartiallyOrderedSet{Polymake.BasicDecoration}(;
    ADJACENCY=pm_object(g),
    DECORATION=dec,
    TOP_NODE=length(sortednodes) - 1,
    BOTTOM_NODE=0,
  )
  opos = PartiallyOrderedSet(pos)
  opos.artificial_bottom = abottom
  opos.artificial_top = atop
  return opos
end

@doc raw"""
    comparability_graph(pos::PartiallyOrderedSet)

Return an undirected graph which has an edge for each pair of comparable
elements in a partially ordered set.

# Examples
```jldoctest
julia> p = face_poset(simplex(2))
Partially ordered set of rank 3 on 8 elements

julia> comparability_graph(p)
Undirected graph with 8 nodes and the following edges:
(2, 1)(3, 1)(4, 1)(5, 1)(5, 2)(5, 4)(6, 1)(6, 2)(6, 3)(7, 1)(7, 3)(7, 4)(8, 1)(8, 2)(8, 3)(8, 4)(8, 5)(8, 6)(8, 7)
```
"""
comparability_graph(pos::PartiallyOrderedSet) =
  Graph{Undirected}(pm_object(pos).COMPARABILITY_GRAPH::Polymake.Graph{Undirected})

# this will compute the comparability graph on first use
# the graph will be kept in the polymake object
function compare(pos::PartiallyOrderedSet, a::Int, b::Int)
  cg = comparability_graph(pos)
  has_edge(cg, a, b) || return false
  return rank(pos, a) < rank(pos, b)
end

function _get_decoration(p::PartiallyOrderedSet, i::Int)
  dec = pm_object(p).DECORATION::Polymake.NodeMap{Directed,Polymake.BasicDecoration}
  elemdec = Polymake._get_entry(dec, i - 1)::Polymake.BasicDecoration
  return Polymake.decoration_face(elemdec)
end

function Base.:<(a::PartiallyOrderedSetElement, b::PartiallyOrderedSetElement)
  @req parent(a) === parent(b) "cannot compare elements in different posets"
  return compare(parent(a), node_id(a), node_id(b))
end

function Base.:(==)(a::PartiallyOrderedSetElement, b::PartiallyOrderedSetElement)
  @req parent(a) === parent(b) "cannot compare elements in different posets"
  return node_id(a) == node_id(b)
end

function Base.hash(pe::PartiallyOrderedSetElement, h::UInt)
  return hash(node_id(pe), hash(parent(pe), h))
end

parent(pose::PartiallyOrderedSetElement) = pose.parent
node_id(pose::PartiallyOrderedSetElement) = pose.node_id

function data(pose::PartiallyOrderedSetElement)
  dec = _get_decoration(parent(pose), node_id(pose))
  # if the poset has the default polymake subsets as node data
  # then we need to shift the indices
  if isdefined(parent(pose), :atomlabels)
    dec = Polymake.to_one_based_indexing(dec)
  end
  return Vector{Int}(dec)
end

Base.length(p::PartiallyOrderedSet) =
  pm_object(p).N_NODES::Int - p.artificial_top - p.artificial_bottom

@doc raw"""
    rank(pe::PartiallyOrderedSetElement)

Return the rank of an element in a partially ordered set.

# Examples
```jldoctest
julia> p = face_poset(simplex(2))
Partially ordered set of rank 3 on 8 elements

julia> rank(greatest_element(p))
3

julia> rank(first(coatoms(p)))
2
```
"""
function rank(pe::PartiallyOrderedSetElement)
  return rank(parent(pe), node_id(pe))
end

function rank(p::PartiallyOrderedSet, i::Int64)
  dec = pm_object(p).DECORATION::Polymake.NodeMap{Directed,Polymake.BasicDecoration}
  elemdec = Polymake._get_entry(dec, i - 1)::Polymake.BasicDecoration
  return Polymake.decoration_rank(elemdec)::Int
end

@doc raw"""
    rank(p::PartiallyOrderedSet)

Return the rank of a partially ordered set, i.e., the rank of a maximal element.

# Examples
```jldoctest
julia> p = face_poset(simplex(2))
Partially ordered set of rank 3 on 8 elements

julia> rank(p)
3
```
"""
function rank(p::PartiallyOrderedSet)
  return rank(p, _top_node(p)) - p.artificial_top
end

@doc raw"""
    maximal_chains([::Type{Int},] p::PartiallyOrderedSet)

Return the maximal chains of a partially ordered set as a vector of sets of elements (or node ids of `Int` is passed as first argument).

# Examples
```jldoctest
julia> pos = face_poset(cube(2))
Partially ordered set of rank 3 on 10 elements

julia> maximal_chains(Int, pos)
8-element Vector{Vector{Int64}}:
 [1, 2, 6, 10]
 [1, 2, 8, 10]
 [1, 3, 7, 10]
 [1, 3, 8, 10]
 [1, 4, 6, 10]
 [1, 4, 9, 10]
 [1, 5, 7, 10]
 [1, 5, 9, 10]

julia> pos = face_poset(simplex(2))
Partially ordered set of rank 3 on 8 elements

julia> maximal_chains(pos)
6-element Vector{Vector{Oscar.PartiallyOrderedSetElement}}:
 [Poset element <Int64[]>, Poset element <[1]>, Poset element <[1, 3]>, Poset element <[1, 2, 3]>]
 [Poset element <Int64[]>, Poset element <[1]>, Poset element <[1, 2]>, Poset element <[1, 2, 3]>]
 [Poset element <Int64[]>, Poset element <[2]>, Poset element <[1, 2]>, Poset element <[1, 2, 3]>]
 [Poset element <Int64[]>, Poset element <[2]>, Poset element <[2, 3]>, Poset element <[1, 2, 3]>]
 [Poset element <Int64[]>, Poset element <[3]>, Poset element <[1, 3]>, Poset element <[1, 2, 3]>]
 [Poset element <Int64[]>, Poset element <[3]>, Poset element <[2, 3]>, Poset element <[1, 2, 3]>]

```
"""
function maximal_chains(p::PartiallyOrderedSet)
  ids = maximal_chains(Int, p)
  return map(x -> Oscar.PartiallyOrderedSetElement.(Ref(p), x), ids)
end

function maximal_chains(::Type{Int}, p::PartiallyOrderedSet)
  mc = Polymake.graph.maximal_chains_of_lattice(
    pm_object(p); ignore_top_node=p.artificial_top, ignore_bottom_node=p.artificial_bottom
  )::IncidenceMatrix
  return Vector.(Polymake.to_one_based_indexing.(Polymake._row.(Ref(mc), 1:nrows(mc))))
end

function Base.show(io::IO, p::PartiallyOrderedSet)
  if is_terse(io)
    print(io, "Partially ordered set")
  else
    print(io, "Partially ordered set of rank $(rank(p)) on $(length(p)) elements")
  end
end

function Base.show(io::IO, pe::PartiallyOrderedSetElement)
  if is_terse(io)
    print(io, "<$(data(pe))>")
  else
    print(io, "Poset element <$(data(pe))>")
  end
end

@doc raw"""
    visualize(p::PartiallyOrderedSet; AtomLabels=[], filename=nothing)

Visualize a partially ordered set.
The labels will be either the ids from the original input data or integer sets
when the poset was created from inclusion relations.
For inclusion based partially ordered sets the keyword `AtomLabels` can be used to specify
a vector of strings with one label per atom to override the default `1:n_atoms(p)` labeling.
The labels of the least and greatest elements are not shown.
The `filename` keyword argument allows writing TikZ visualization code to `filename`.

!!! note
    This will always show the greatest and least elements even if it does not correspond to an element of the partially ordered set.
"""
function visualize(
  p::PartiallyOrderedSet; filename::Union{Nothing,String}=nothing, kwargs...
)
  if isdefined(p, :atomlabels) && !haskey(kwargs, :AtomLabels)
    kwargs = (kwargs..., :AtomLabels => p.atomlabels)
  end
  vpos = Polymake.visual(Polymake.Visual, pm_object(p); kwargs...)
  # todo: check if we can use svg in jupyter notebooks
  if !isnothing(filename)
    Polymake.call_function(Nothing, :graph, :tikz, vpos; File=filename)
  else
    Polymake.call_function(Nothing, :graph, :tikz, vpos;)
  end
end

@doc raw"""
    face_poset(p::Union{Polyhedron,Cone,PolyhedralFan,PolyhedralComplex,SimplicialComplex})

Return a partially ordered set encoding the face inclusion relations of a polyhedral object.

Note that a polyhedron has the entire polyhedron as the greatest element while the normal fan does not have a greatest element.

# Examples
```jldoctest
julia> p = face_poset(dodecahedron())
Partially ordered set of rank 4 on 64 elements

julia> pf = face_poset(normal_fan(dodecahedron()))
Partially ordered set of rank 3 on 63 elements
```
"""
function face_poset(
  p::Union{Polyhedron,Cone,PolyhedralFan,PolyhedralComplex,SimplicialComplex}
)
  pos = partially_ordered_set(pm_object(p).HASSE_DIAGRAM)
  pos.atomlabels = collect(1:n_atoms(pos))
  pos.artificial_top = p isa Union{PolyhedralFan,PolyhedralComplex,SimplicialComplex}
  return pos
end

@doc raw"""
    lattice_of_cyclic_flats(m::Matroid)

Return the lattice of cyclic flats of a matroid.

# Examples
```jldoctest
julia> p = lattice_of_cyclic_flats(fano_matroid())
Partially ordered set of rank 3 on 9 elements
```
"""
function lattice_of_cyclic_flats(m::Matroid)
  pos = partially_ordered_set(pm_object(m).LATTICE_OF_CYCLIC_FLATS)
  pos.atomlabels = string.(matroid_groundset(m))
  return pos
end

@doc raw"""
    lattice_of_flats(m::Matroid)

Return the lattice of flats of a matroid.

# Examples
```jldoctest
julia> p = lattice_of_flats(fano_matroid())
Partially ordered set of rank 3 on 16 elements
```
"""
function lattice_of_flats(m::Matroid)
  pos = partially_ordered_set(pm_object(m).LATTICE_OF_FLATS)
  pos.atomlabels = string.(matroid_groundset(m))
  return pos
end

@doc raw"""
    order_polytope(p::PartiallyOrderedSet)

Return the order polytope of a poset, see [Sta86](@cite).

# Examples
```jldoctest
julia> p = face_poset(simplex(2))
Partially ordered set of rank 3 on 8 elements

julia> op = order_polytope(p)
Polytope in ambient dimension 6

julia> f_vector(op)
6-element Vector{ZZRingElem}:
 18
 73
 129
 116
 54
 12
```
"""
function order_polytope(p::PartiallyOrderedSet)
  return polyhedron(Polymake.polytope.order_polytope(pm_object(p)))
end

@doc raw"""
    chain_polytope(p::PartiallyOrderedSet)

Return the chain polytope of a poset, see [Sta86](@cite).

# Examples
```jldoctest
julia> p = face_poset(cube(2))
Partially ordered set of rank 3 on 10 elements

julia> cp = chain_polytope(p)
Polytope in ambient dimension 8
```
"""
function chain_polytope(p::PartiallyOrderedSet)
  return polyhedron(Polymake.polytope.chain_polytope(pm_object(p)))
end

@doc raw"""
    maximal_ranked_poset(v::AbstractVector{Int})

Return a maximal ranked partially ordered set with number of nodes per rank given in
`v`, from bottom to top and excluding the least and greatest elements.

See [AFJ25](@cite) for details.

# Examples
```jldoctest
julia> pos = maximal_ranked_poset([2,4,3])
Partially ordered set of rank 4 on 11 elements
```
"""
function maximal_ranked_poset(v::AbstractVector{Int})
  pos = partially_ordered_set(Polymake.graph.maximal_ranked_poset(Polymake.Array{Int}(v)))
  pos.artificial_top = false
  pos.artificial_bottom = false
  return pos
end

@doc raw"""
    graph(p::PartiallyOrderedSet)

Return the directed graph of covering relations in a partially ordered set.
To keep the node ids consistent this will include artificial top and bottom nodes.

# Examples
```jldoctest
julia> pos = maximal_ranked_poset([2,4,3])
Partially ordered set of rank 4 on 11 elements

julia> graph(pos)
Directed graph with 11 nodes and the following edges:
(1, 2)(1, 3)(2, 4)(2, 5)(2, 6)(2, 7)(3, 4)(3, 5)(3, 6)(3, 7)(4, 8)(4, 9)(4, 10)(5, 8)(5, 9)(5, 10)(6, 8)(6, 9)(6, 10)(7, 8)(7, 9)(7, 10)(8, 11)(9, 11)(10, 11)
```
"""
graph(p::PartiallyOrderedSet) = Graph{Directed}(pm_object(p).ADJACENCY)

@doc raw"""
    node_ranks(p::PartiallyOrderedSet)

Return a dictionary mapping each node index to its rank.
"""
node_ranks(p::PartiallyOrderedSet) = Dict((i => rank(p, i)) for i in 1:length(p))

@doc raw"""
    n_atoms(p::PartiallyOrderedSet)

Return the number of atoms in a partially ordered set.

# Examples
```jldoctest
julia> pos = maximal_ranked_poset([2,4,3])
Partially ordered set of rank 4 on 11 elements

julia> n_atoms(pos)
2
```
"""
n_atoms(p::PartiallyOrderedSet) = length(atoms(p))

@doc raw"""
    n_coatoms(p::PartiallyOrderedSet)

Return the number of co-atoms in a partially ordered set.

# Examples
```jldoctest
julia> pos = maximal_ranked_poset([2,4,3])
Partially ordered set of rank 4 on 11 elements

julia> n_coatoms(pos)
3
```
"""
n_coatoms(p::PartiallyOrderedSet) = length(coatoms(p))

@doc raw"""
    atoms(p::PartiallyOrderedSet)

Return the atoms in a partially ordered set.

# Examples
```jldoctest
julia> pos = maximal_ranked_poset([2,4,3])
Partially ordered set of rank 4 on 11 elements

julia> atoms(pos)
2-element Vector{Oscar.PartiallyOrderedSetElement}:
 Poset element <[1]>
 Poset element <[2]>
```
"""
atoms(p::PartiallyOrderedSet) =
  PartiallyOrderedSetElement.(Ref(p), outneighbors(graph(p), node_id(least_element(p))))

@doc raw"""
    coatoms(p::PartiallyOrderedSet)

Return the co-atoms in a partially ordered set.

# Examples
```jldoctest
julia> pos = maximal_ranked_poset([2,4,3])
Partially ordered set of rank 4 on 11 elements

julia> coatoms(pos)
3-element Vector{Oscar.PartiallyOrderedSetElement}:
 Poset element <[7]>
 Poset element <[8]>
 Poset element <[9]>
```
"""
coatoms(p::PartiallyOrderedSet) =
  PartiallyOrderedSetElement.(Ref(p), inneighbors(graph(p), node_id(greatest_element(p))))

_bottom_node(p::PartiallyOrderedSet) = 1 + pm_object(p).BOTTOM_NODE::Int
_top_node(p::PartiallyOrderedSet) = 1 + pm_object(p).TOP_NODE::Int

function _ids_of_rank(p::PartiallyOrderedSet, rk::Int)
  pv = Polymake.graph.nodes_of_rank(pm_object(p), rk)
  ids = Polymake.@convert_to Array{Int} pv
  return Polymake.to_one_based_indexing(ids)
end

@doc raw"""
    least_element(p::PartiallyOrderedSet)

Return the unique minimal element in a partially ordered set, if it exists.
Throws an error otherwise.

# Examples
```jldoctest
julia> pos = maximal_ranked_poset([2,4,3])
Partially ordered set of rank 4 on 11 elements

julia> least_element(pos)
Poset element <[0]>
```
"""
function least_element(p::PartiallyOrderedSet)
  bottom = _bottom_node(p)
  if p.artificial_bottom
    br = rank(p, bottom)
    nb = _ids_of_rank(p, br + 1)
    @req length(nb) == 1 "No unique minimal element"
    bottom = first(nb)
  end
  return PartiallyOrderedSetElement(p, bottom)
end

@doc raw"""
    greatest_element(p::PartiallyOrderedSet)

Return the unique maximal element in a partially ordered set, if it exists.
Throws an error otherwise.

# Examples
```jldoctest
julia> pos = maximal_ranked_poset([2,4,3])
Partially ordered set of rank 4 on 11 elements

julia>  greatest_element(pos)
Poset element <[10]>
```
"""
function greatest_element(p::PartiallyOrderedSet)
  top = _top_node(p)
  if p.artificial_top
    tr = rank(p, top)
    nb = _ids_of_rank(p, tr - 1)
    @req length(nb) == 1 "No unique maximal element"
    top = first(nb)
  end
  return PartiallyOrderedSetElement(p, top)
end

@doc raw"""
    element(p::PartiallyOrderedSet, i::Int)

Return the element in a partially ordered set with give node id.

# Examples
```jldoctest
julia> fl = face_poset(cube(2))
Partially ordered set of rank 3 on 10 elements

julia> element(fl,7)
Poset element <[2, 4]>
```
"""
function element(p::PartiallyOrderedSet, i::Int)
  @req _bottom_node(p) + p.artificial_bottom <= i <= _top_node(p) - p.artificial_top "invalid node id"
  return PartiallyOrderedSetElement(p, i)
end

@doc raw"""
    elements(p::PartiallyOrderedSet)

Return all elements in a partially ordered set.

# Examples
```jldoctest
julia> pos = maximal_ranked_poset([2,4,3])
Partially ordered set of rank 4 on 11 elements

julia> length(elements(pos))
11
```
"""
elements(p::PartiallyOrderedSet) =
  PartiallyOrderedSetElement.(
    Ref(p), (_bottom_node(p) + p.artificial_bottom):(_top_node(p) - p.artificial_top)
  )

@doc raw"""
    elements_of_rank(p::PartiallyOrderedSet, rk::Int)

Return the elements in a partially ordered set with a given rank.

# Examples
```jldoctest
julia> pos = maximal_ranked_poset([2,4,3])
Partially ordered set of rank 4 on 11 elements

julia> length(elements_of_rank(pos, 2))
4
```
"""
elements_of_rank(p::PartiallyOrderedSet, rk::Int) =
  PartiallyOrderedSetElement.(Ref(p), _ids_of_rank(p, rk))
