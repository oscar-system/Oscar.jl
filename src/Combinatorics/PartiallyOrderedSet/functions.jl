function pm_object(P::PartiallyOrderedSet)
    return P.pm_poset
end

function partially_ordered_set(pp::Polymake.BigObject)
  @req startswith(Polymake.type_name(pp), "PartiallyOrderedSet") "Not a polymake PartiallyOrderedSet: $(Polymake.type_name(pp))"
  return PartiallyOrderedSet(pp)
end

@doc raw"""
    partially_ordered_set(covrels::Matrix{Int})
￼
Construct a partially ordered set from covering relations `covrels`, given in
the form of the adjacency matrix of a Hasse diagram. The covering relations
must be given in topological order, i.e. `covrels` must be strictly upper
triangular.
"""
function partially_ordered_set(covrels::Matrix{Int})
  pg = Polymake.Graph{Directed}(IncidenceMatrix(covrels))
  nelem = nrows(covrels)
  return partially_ordered_set(Graph{Directed}(pg), 1, nelem)
end


@doc raw"""
￼   partially_ordered_set_from_inclusions(i::IncidenceMatrix)
￼
Construct an inclusion based partially ordered with rows of `i` as co-atoms.
"""
function partially_ordered_set_from_inclusions(i::IncidenceMatrix)
  # we want to make sure it is always built from the bottom up
  # and just need to specify some arbitrary upper bound for the rank
  pos = Polymake.call_function(:polytope, :lower_hasse_diagram, i, ncols(i))
  return partially_ordered_set(pos)
end


_pmdec(elem::Int, r::Int) = Polymake.BasicDecoration(Polymake.Set{Int}(elem), r)
_pmdec(r::Int) = Polymake.BasicDecoration(Polymake.Set{Int}(), r)

@doc raw"""
￼   partially_ordered_set(g::Graph{Directed}, minimal_element::Int, maximal_element::Int)
￼
Construct a partially ordered set from a directed graph describing the Hasse diagram.
The graph must be acyclic and have unique minimal and maximal elements. Edges must be
directed from bottom to top.
"""
function partially_ordered_set(g::Graph{Directed}, minimal_element::Int, maximal_element::Int)
  stack = [minimal_element]
  gc = Polymake.Graph{Directed}(pm_object(g))
  dec = Polymake.NodeMap{Directed, Polymake.BasicDecoration}(gc)
  # we need to generate some valid rank function
  es = Polymake.Set{Int}()
  Polymake._set_entry(dec, minimal_element-1, _pmdec(minimal_element, 0))
  pos = 1
  # we keep the elements in the stack
  while pos <= length(stack)
    src = stack[pos]
    r = Polymake.decoration_rank(Polymake._get_entry(dec, src-1))
    for t in neighbors(g, src)
      rr = Polymake.decoration_rank(Polymake._get_entry(dec, t-1))
      Polymake._set_entry(dec, t-1, _pmdec(t, max(r+1, rr)))
      push!(stack, t)
    end
    pos += 1
  end
  println(stack)
  println(gc)
  Polymake.common.permute_graph(gc, Polymake.to_zero_based_indexing(stack))
  pos = Polymake.graph.PartiallyOrderedSet{Polymake.BasicDecoration}(
                                                                     ADJACENCY=gc,
                                                                     DECORATION=dec,
                                                                     TOP_NODE=length(stack)-1,
                                                                     BOTTOM_NODE=0
                                                                    )
  return PartiallyOrderedSet(pos)
end

@doc raw"""
    partially_ordered_set(g::Graph{Directed}, node_ranks::Dict{Int,Int})
￼
Construct a partially ordered set from a directed graph describing the Hasse diagram.
The graph must be acyclic and have unique minimal and maximal elements.
A The dictionary `node_ranks` must give a valid rank for each node in the graph, strictly
increasing from the unique minimal element to the unique maximal element.
"""
function partially_ordered_set(g::Graph{Directed}, node_ranks::Dict{Int,Int})
  gc = Polymake.Graph{Directed}(pm_object(g))
  dec = Polymake.NodeMap{Directed, Polymake.BasicDecoration}(gc)
  for n in 1:n_vertices(g)
    Polymake._set_entry(dec, n-1, _pmdec(n, node_ranks[n]))
  end
  bottom = findmin(node_ranks)[2]
  top = findmax(node_ranks)[2]
  pos = Polymake.graph.PartiallyOrderedSet{Polymake.BasicDecoration}(
                                                                     ADJACENCY=gc,
                                                                     DECORATION=dec,
                                                                     TOP_NODE=top-1,
                                                                     BOTTOM_NODE=bottom-1
                                                                    )
  return PartiallyOrderedSet(pos)
end

@doc raw"""
    comparability_graph(pos::PartiallyOrderedSet)

Return an undirected graph which has an edge for each pair of comparable
elements in a partially ordered set.
"""
comparability_graph(pos::PartiallyOrderedSet) = Graph{Undirected}(pm_object(pos).COMPARABILITY_GRAPH::Polymake.Graph{Undirected})

# this will compute the comparability graph on first use
# the graph will be kept in the polymake object
function compare(pos::PartiallyOrderedSet, a::Int, b::Int)
  cg = comparability_graph(pos)
  has_edge(cg, a, b) || return false
  return rank(pos, a) < rank(pos, b)
end

function compare(a::PartiallyOrderedSetElement, b::PartiallyOrderedSetElement)
  @req parent(a) === parent(b) "cannot compare elements in different posets"
  return compare(parent(a), data(a), data(b))
end

parent(pose::PartiallyOrderedSetElement) = pose.parent
data(pose::PartiallyOrderedSetElement) = pose.node_id

Base.length(p::PartiallyOrderedSet) = pm_object(p).N_NODES::Int

@doc raw"""
    rank(p::PartiallyOrderedSet, i::Int64)

Return the rank of an element with index `i` in a partially ordered set.
"""
function rank(p::PartiallyOrderedSet, i::Int64)
  dec = pm_object(p).DECORATION::Polymake.NodeMap{Directed,Polymake.BasicDecoration}
  elemdec = Polymake._get_entry(dec, i-1)::Polymake.BasicDecoration
  return Polymake.decoration_rank(elemdec)::Int
end

@doc raw"""
    maximal_chains(p::PartiallyOrderedSet)

Return the maximal chains of a partially ordered set as a vector of sets.
The nodes in the chains are not ordered.
"""
function maximal_chains(p::PartiallyOrderedSet)
  mc = Polymake.graph.maximal_chains_of_lattice(pm_object(p))::IncidenceMatrix
  return row.(Ref(mc), 1:nrows(mc))
end

function Base.show(io::IO, p::PartiallyOrderedSet)
  if is_terse(io)
    print(io, "Partially ordered set")
  else
    print(io, "Partially ordered set on $(length(p)) elements")
  end
end

function visualize(p::PartiallyOrderedSet; kwargs...)
  Polymake.visual(pm_object(p); kwargs...)
end

@doc raw"""
    face_lattice(p::Union{Polyhedron,Cone,PolyhedralFan,PolyhedralComplex})

Return a partially ordered set encoding the face inclusion relations of a polyhedral object.
"""
function face_lattice(p::Union{Polyhedron,Cone,PolyhedralFan,PolyhedralComplex})
  return partially_ordered_set(pm_object(p).HASSE_DIAGRAM)
end

@doc raw"""
    order_polytope(p::PartiallyOrderedSet)

Return the order polytope of a poset.
See Stanley, Discr Comput Geom 1 (1986).
TODO: ref
"""
function order_polytope(p::PartiallyOrderedSet)
  return polyhedron(Polymake.polytope.order_polytope(pm_object(p)))
end

@doc raw"""
    chain_polytope(p::PartiallyOrderedSet)

Return the chain polytope of a poset.
See Stanley, Discr Comput Geom 1 (1986).
TODO: ref
"""
function chain_polytope(p::PartiallyOrderedSet)
  return polyhedron(Polymake.polytope.chain_polytope(pm_object(p)))
end

@doc raw"""
    maximal_ranked_poset(v::AbstractVector{Int})

Maximal ranked partially ordered set with number of nodes per rank given in
`v`, from bottom to top and excluding the minimal and maximal elements.

See Ahmad, Fourier, Joswig, arXiv:2309.01626
TODO: ref
"""
function maximal_ranked_poset(v::AbstractVector{Int})
  return partially_ordered_set(Polymake.graph.maximal_ranked_poset(Polymake.Array{Int}(v)))
end

@doc raw"""
    graph(p::PartiallyOrderedSet)

Return the directed graph of covering relations in a partially ordered set.
"""
graph(p::PartiallyOrderedSet) = Graph{Directed}(pm_object(p).ADJACENCY)

@doc raw"""
    node_ranks(p::PartiallyOrderedSet)

Return a dictionary mapping each node index to its rank.
"""
node_ranks(p::PartiallyOrderedSet) = Dict((i => rank(p, i)) for i in 1:length(p))
