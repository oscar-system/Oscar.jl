function pm_object(P::PartiallyOrderedSet)
    return P.pm_poset
end

function partially_ordered_set(pp::Polymake.BigObject)
  @req startswith(Polymake.type_name(pp), "PartiallyOrderedSet") "Not a polymake PartiallyOrderedSet: $(Polymake.type_name(pp))"
  return PartiallyOrderedSet(pp)
end

@doc raw"""
    partially_ordered_set(covrels::Matrix{Int})

Construct a partially ordered set from covering relations `covrels`, given in
the form of the adjacency matrix of a Hasse diagram. The covering relations
must be given in topological order, i.e., `covrels` must be strictly upper
triangular.

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

The construction only distinguishes zero coefficients (to indicate the absence) and nonzero coefficients (to indicate the presence of a covering relation).
"""
function partially_ordered_set(covrels::Matrix{Int})
  pg = Polymake.Graph{Directed}(IncidenceMatrix(covrels))
  nelem = nrows(covrels)
  return partially_ordered_set(Graph{Directed}(pg))
end


@doc raw"""
    partially_ordered_set_from_inclusions(I::IncidenceMatrix)

Construct an inclusion based partially ordered with rows of `I` as co-atoms.
That is, the elements of the poset are the given sets and their intersections.

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
The graph must be acyclic and have unique minimal and maximal elements. Edges must be
directed from minimal to maximal element.

# Examples
```jldoctest
julia> g = graph_from_edges(Directed, [[2,1],[2,4],[1,3],[4,3]])
Directed graph with 4 nodes and the following edges:
(1, 3)(2, 1)(2, 4)(4, 3)

julia> pos = partially_ordered_set(g)
Partially ordered set of rank 2 on 4 elements
```
"""
function partially_ordered_set(g::Graph{Directed})
  @req n_vertices(g) >= 2  "Graph must have at least two nodes"
  gc = Polymake.Graph{Directed}(pm_object(g))
  dec = Polymake.NodeMap{Directed, Polymake.BasicDecoration}(gc)
  _, least_element = findmax(degree(g))
  deg, _ = findmin(degree(g))
  @req deg == 0 "no source element found"
  queue = [least_element]
  Polymake._set_entry(dec, least_element-1, _pmdec(least_element, 0))
  pos = 1
  # we need to generate a valid rank function
  visited = BitSet(queue)
  # we keep the elements in the queue for permuting the graph at the end
  while pos <= length(queue)
    src = queue[pos]
    r = Polymake.decoration_rank(Polymake._get_entry(dec, src-1))
    for t in neighbors(g, src)
      rr = Polymake.decoration_rank(Polymake._get_entry(dec, t-1))
      Polymake._set_entry(dec, t-1, _pmdec(t, max(r+1, rr)))
      if !(t in visited)
        push!(queue, t)
        push!(visited, t)
      end
    end
    pos += 1
  end
  Polymake.call_function(:common, :permute_graph, gc, Polymake.to_zero_based_indexing(queue))
  pos = Polymake.graph.PartiallyOrderedSet{Polymake.BasicDecoration}(
                                                                     ADJACENCY=gc,
                                                                     DECORATION=dec,
                                                                     TOP_NODE=length(queue)-1,
                                                                     BOTTOM_NODE=0
                                                                    )
  return PartiallyOrderedSet(pos)
end

@doc raw"""
    partially_ordered_set(g::Graph{Directed}, node_ranks::Dict{Int,Int})

Construct a partially ordered set from a directed graph describing the Hasse diagram.
The graph must be acyclic and have unique minimal and maximal elements.
The dictionary `node_ranks` must give a valid rank for each node in the graph, strictly
increasing from the unique minimal element to the unique maximal element.
The rank difference between two adjacent nodes may be larger than one.

# Examples
```jldoctest
julia> g = graph_from_edges(Directed, [[2,1],[2,4],[1,3],[4,3],[3,6],[2,5],[5,6]])
Directed graph with 6 nodes and the following edges:
(1, 3)(2, 1)(2, 4)(2, 5)(3, 6)(4, 3)(5, 6)

julia> pos = partially_ordered_set(g, Dict(2=>0, 1=>1, 4=>2, 3=>3, 5=>2, 6=>4))
Partially ordered set of rank 4 on 6 elements
```
"""
function partially_ordered_set(g::Graph{Directed}, node_ranks::Dict{Int,Int})
  @req n_vertices(g) >= 2  "Graph must have at least two nodes"
  gc = Polymake.Graph{Directed}(pm_object(g))
  dec = Polymake.NodeMap{Directed, Polymake.BasicDecoration}(gc)
  for n in 1:n_vertices(g)
    Polymake._set_entry(dec, n-1, _pmdec(n, node_ranks[n]))
  end
  sorted = first.(sort(collect(node_ranks), by=last))
  Polymake.call_function(:common, :permute_graph, gc, Polymake.to_zero_based_indexing(sorted))
  pos = Polymake.graph.PartiallyOrderedSet{Polymake.BasicDecoration}(
                                                                     ADJACENCY=gc,
                                                                     DECORATION=dec,
                                                                     TOP_NODE=length(sorted)-1,
                                                                     BOTTOM_NODE=0
                                                                    )
  return PartiallyOrderedSet(pos)
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
comparability_graph(pos::PartiallyOrderedSet) = Graph{Undirected}(pm_object(pos).COMPARABILITY_GRAPH::Polymake.Graph{Undirected})

# this will compute the comparability graph on first use
# the graph will be kept in the polymake object
function compare(pos::PartiallyOrderedSet, a::Int, b::Int)
  cg = comparability_graph(pos)
  has_edge(cg, a, b) || return false
  return rank(pos, a) < rank(pos, b)
end

function _get_decoration(p::PartiallyOrderedSet, i::Int)
  dec = pm_object(p).DECORATION::Polymake.NodeMap{Directed,Polymake.BasicDecoration}
  elemdec = Polymake._get_entry(dec, i-1)::Polymake.BasicDecoration
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

Base.length(p::PartiallyOrderedSet) = pm_object(p).N_NODES::Int - p.artificial_top

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
  elemdec = Polymake._get_entry(dec, i-1)::Polymake.BasicDecoration
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
    maximal_chains(p::PartiallyOrderedSet)

Return the maximal chains of a partially ordered set as a vector of sets of node ids.

# Examples
```jldoctest
julia> pos = face_poset(cube(2))
Partially ordered set of rank 3 on 10 elements

julia> maximal_chains(pos)
8-element Vector{Vector{Int64}}:
 [1, 2, 6, 10]
 [1, 2, 8, 10]
 [1, 3, 7, 10]
 [1, 3, 8, 10]
 [1, 4, 6, 10]
 [1, 4, 9, 10]
 [1, 5, 7, 10]
 [1, 5, 9, 10]
```
"""
function maximal_chains(p::PartiallyOrderedSet)
  mc = Polymake.graph.maximal_chains_of_lattice(pm_object(p))::IncidenceMatrix
  # TODO check with top and bottom
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
    This will always show the greatest element even if it does not correspond to an element of the partially ordered set.
"""
function visualize(p::PartiallyOrderedSet; filename::Union{Nothing, String}=nothing, kwargs...)
  if isdefined(p, :atomlabels) && !haskey(kwargs, :AtomLabels)
    kwargs = (kwargs..., :AtomLabels=>p.atomlabels)
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
function face_poset(p::Union{Polyhedron,Cone,PolyhedralFan,PolyhedralComplex,SimplicialComplex})
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

Maximal ranked partially ordered set with number of nodes per rank given in
`v`, from bottom to top and excluding the least and greatest elements.

See [AFJ25](@cite) for details.

# Examples
```jldoctest
julia> pos = maximal_ranked_poset([2,4,3])
Partially ordered set of rank 4 on 11 elements
```
"""
function maximal_ranked_poset(v::AbstractVector{Int})
  return partially_ordered_set(Polymake.graph.maximal_ranked_poset(Polymake.Array{Int}(v)))
end

@doc raw"""
    graph(p::PartiallyOrderedSet)

Return the directed graph of covering relations in a partially ordered set.

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

# Examples
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
atoms(p::PartiallyOrderedSet) = PartiallyOrderedSetElement.(Ref(p), outneighbors(graph(p), node_id(least_element(p))))

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
coatoms(p::PartiallyOrderedSet) = PartiallyOrderedSetElement.(Ref(p), inneighbors(graph(p), node_id(greatest_element(p))))

_bottom_node(p::PartiallyOrderedSet) = 1+pm_object(p).BOTTOM_NODE::Int
_top_node(p::PartiallyOrderedSet) = 1+pm_object(p).TOP_NODE::Int

function _ids_of_rank(p::PartiallyOrderedSet, rk::Int)
  pv = Polymake.graph.nodes_of_rank(pm_object(p), rk)
  ids = Polymake.@convert_to Array{Int} pv
  return Polymake.to_one_based_indexing(ids)
end

@doc raw"""
    least_element(p::PartiallyOrderedSet)

Return the unique minimal element in a partially ordered set.

# Examples
```jldoctest
julia> pos = maximal_ranked_poset([2,4,3])
Partially ordered set of rank 4 on 11 elements

julia> least_element(pos)
Poset element <[0]>
```
"""
least_element(p::PartiallyOrderedSet) = PartiallyOrderedSetElement(p, _bottom_node(p))

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
    nb = _ids_of_rank(p, tr-1)
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
element(p::PartiallyOrderedSet, i::Int) = PartiallyOrderedSetElement(p, i)

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
elements(p::PartiallyOrderedSet) = PartiallyOrderedSetElement.(Ref(p), 0:n_vertices(graph(p)) - 1 - p.artificial_top)

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
elements_of_rank(p::PartiallyOrderedSet, rk::Int) = PartiallyOrderedSetElement.(Ref(p), _ids_of_rank(p, rk))
