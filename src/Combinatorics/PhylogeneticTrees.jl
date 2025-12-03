struct PhylogeneticTree{T <: Union{Float64, QQFieldElem}} <: AbstractGraph{Directed}
  pm_ptree::Polymake.LibPolymake.BigObjectAllocated
  vertex_perm::Vector{Int}
end

function PhylogeneticTree{T}(pm::Polymake.LibPolymake.BigObjectAllocated) where T <: Union{Float64, QQFieldElem}
  return PhylogeneticTree{T}(pm, collect(1:pm.N_NODES))
end


function pm_object(PT::PhylogeneticTree)
  return PT.pm_ptree
end

vertex_perm(pt::PhylogeneticTree) = pt.vertex_perm

@doc raw"""
    phylogenetic_tree(T::Type{<:Union{Float64, QQFieldElem}}, newick::String)
    phylogenetic_tree(newick::String)
Construct a rooted phylogenetic tree with Newick representation `newick`.
`T` indicates the numerical type of the edge lengths.

Calling `phylogenetic_tree` without `T` will default to using QQFieldElem.

# Examples
Make a phylogenetic tree with 4 leaves from its Newick representation and print
its taxa and cophenetic matrix.
```jldoctest
julia> phylo_t = phylogenetic_tree(Float64, "((H:3,(C:1,B:1):2):1,G:4);");

julia> taxa(phylo_t)
4-element Vector{String}:
 "B"
 "C"
 "G"
 "H"

julia> cophenetic_matrix(phylo_t)
4×4 Matrix{Float64}:
 0.0  2.0  8.0  6.0
 2.0  0.0  8.0  6.0
 8.0  8.0  0.0  8.0
 6.0  6.0  8.0  0.0

julia> typeof(phylogenetic_tree("((H:3,(C:1,B:1):2):1,G:4);"))
Oscar.PhylogeneticTree{QQFieldElem}
```
"""
function phylogenetic_tree(T::Type{<:Union{Float64, QQFieldElem}}, newick::String)
  pm_ptree = Polymake.graph.PhylogeneticTree{Polymake.convert_to_pm_type(T)}(NEWICK = newick)
  # load graph properties
  pm_ptree.ADJACENCY

  return PhylogeneticTree{T}(pm_ptree)
end

phylogenetic_tree(newick::String) = phylogenetic_tree(QQFieldElem, newick)

@doc raw"""
    phylogenetic_tree(M::Matrix{T}, taxa::Vector{String}) where T <: Union{Float64, QQFieldElem}

Construct a phylogenetic tree with cophenetic matrix `M` and taxa `taxa`. The matrix `M` must be
ultrametric, otherwise an error will be thrown.

# Examples
Make a phylogenetic tree on 4 taxa with given cophenetic matrix and print one Newick representation.

```jldoctest
julia> mat = [0. 2 8 6; 2 0 8 6; 8 8 0 8; 6 6 8 0]
4×4 Matrix{Float64}:
 0.0  2.0  8.0  6.0
 2.0  0.0  8.0  6.0
 8.0  8.0  0.0  8.0
 6.0  6.0  8.0  0.0

julia> tax = ["Bonobo", "Chimpanzee", "Gorilla", "Human"]
4-element Vector{String}:
 "Bonobo"
 "Chimpanzee"
 "Gorilla"
 "Human"

julia> tree_mat = phylogenetic_tree(mat, tax);

julia> newick(tree_mat)
"Gorilla:4,(Human:3,(Bonobo:1,Chimpanzee:1):2):1;"
```
"""
function phylogenetic_tree(M::Matrix{Float64}, taxa::Vector{String})
  n_taxa = length(taxa)
  @req (n_taxa, n_taxa) == size(M) "Number of taxa should match the rows and columns of the given matrix"
  pm_ptree = Polymake.graph.PhylogeneticTree{Float64}(COPHENETIC_MATRIX = M, TAXA = taxa)
  return PhylogeneticTree{Float64}(pm_ptree)
end

function phylogenetic_tree(M::QQMatrix, taxa::Vector{String})
  n_taxa = length(taxa)
  @req (n_taxa, n_taxa) == size(M) "Number of taxa should match the rows and columns of the given matrix"
  pm_ptree = Polymake.graph.PhylogeneticTree{Rational}(
    COPHENETIC_MATRIX = M, TAXA = taxa
  )
  return PhylogeneticTree{QQFieldElem}(pm_ptree)
end

@doc raw"""
     phylogenetic_tree(T::Type{<:Union{Float64, QQFieldElem}}, g::Graph{Directed}; check=true)
     phylogenetic_tree(g::Graph{Directed}; check=true)

Constructs a phylogenetic tree from a directed graph.
If `check` is set to `true` we check that `g` is a tree. Calling `phylogenetic_tree` without a type `T` will default to `QQFieldelem`.

If the graph is labeled on its leaves using the label name `:leaves` and or if the graph is labeled on its edges with the label name `:distance` then this data will persist to the `PhylogeneticTree`.

# Examples
```jldoctest
julia> tree = graph_from_edges(Directed,[[4,1],[4,2],[4,3], [5, 4], [5, 6]])
Directed graph with 6 nodes and the following edges:
(4, 1)(4, 2)(4, 3)(5, 4)(5, 6)

julia> pt = phylogenetic_tree(tree)
Phylogenetic tree with QQFieldElem type coefficients

julia> newick(pt)
"(leaf 2:1,leaf 3:1,leaf 1:1):1,leaf 4:1;"

julia> label!(tree, Dict((5, 6) => 1.0,
                         (5, 4) => 2.0,
                         (4, 1) => 3.0,
                         (4, 2) => 4.0,
                         (4, 3) => 5.0), nothing; name=:distance)
Directed graph with 6 nodes and the following labeling(s):
label: distance
(4, 1) -> 3.0
(4, 2) -> 4.0
(4, 3) -> 5.0
(5, 4) -> 2.0
(5, 6) -> 1.0

julia> label!(tree, nothing, Dict(1 => "a", 2 => "b", 3 => "c", 6 => "d"); name=:leaves)
Directed graph with 6 nodes and the following labeling(s):
label: leaves
1 -> a
2 -> b
3 -> c
4 ->
5 ->
6 -> d
label: distance
(4, 1) -> 3.0
(4, 2) -> 4.0
(4, 3) -> 5.0
(5, 4) -> 2.0
(5, 6) -> 1.0

julia> pt = phylogenetic_tree(Float64, tree)
Phylogenetic tree with Float64 type coefficients

julia> cophenetic_matrix(pt)
4×4 Matrix{Float64}:
 0.0  7.0  8.0  6.0
 7.0  0.0  9.0  7.0
 8.0  9.0  0.0  8.0
 6.0  7.0  8.0  0.0
```
"""
function phylogenetic_tree(T::Type{<:Union{Float64, QQFieldElem}},
                           g::Graph{Directed};
                           check=true)
  @req !check || _is_tree(g) "Input must be a tree"
  r = root(g)
  p = collect(1:n_vertices(g))
  # root needs to be labeled by 1, se we just transpose 2 vertices
  # for the underlying polymake graph
  p[1], p[r] = r, 1
  undir_g = graph_from_edges(Undirected, edges(g))
  Polymake.call_function(
    :common, :permute_graph, pm_object(undir_g), Polymake.to_zero_based_indexing(p)
  )

  leaves = findall(x -> isone(x[1]) && iszero(x[2]), collect(zip(indegree(g), outdegree(g))))
  lvd = Dict{String, Int}()
  lad = Dict{Int, String}()
  for (i, v) in enumerate(leaves)
    # sending to zero based indexing is important to guarantee cophenetic
    # matrix still works
    if has_attribute(g, :leaves)
      lvd[g.leaves[v]] = Polymake.to_zero_based_indexing(p[v])
      lad[p[v]] = g.leaves[v]
    else
      lvd["leaf $i"] = Polymake.to_zero_based_indexing(p[v])
      lad[p[v]] = "leaf $i"
    end
  end
  lv = Polymake.Map{Polymake.CxxWrap.StdLib.StdString,Int}(lvd)

  eld = Dict{NTuple{2, Int}, T}()
  for e in edges(g)
    s, d = src(e), dst(e)
    if has_attribute(g, :distance)
      eld[p[s], p[d]] = g.distance[s, d]
    else
      eld[p[s], p[d]] = 1
    end
  end

  el = Polymake.EdgeMap{Undirected,Polymake.convert_to_pm_type(T)}(pm_object(undir_g), eld)
  la = Polymake.NodeMap{Undirected,String}(pm_object(undir_g), lad)

  pt = Polymake.graph.PhylogeneticTree{Polymake.convert_to_pm_type(T)}()

  Polymake.take(pt, "ADJACENCY", pm_object(undir_g))
  Polymake._take_graph_map(pt, "ADJACENCY", "EDGE_LENGTHS", el)
  Polymake._take_graph_map(pt, "ADJACENCY", "LABELS", la)
  Polymake.take(pt, "LEAVES", lv)
  Polymake.call_method(pt, :commit)

  return PhylogeneticTree{T}(pt, p)
end

phylogenetic_tree(g::Graph{Directed}; check=true) = phylogenetic_tree(QQFieldElem, g; check=check)

@doc raw"""
    adjacency_tree(ptree::PhylogeneticTree)

Return the underlying graph of the phylogenetic tree `ptree`.

# Examples
Make a phylogenetic tree with given Newick format and print its underlying graph.

```jldoctest
julia> ptree = phylogenetic_tree(Float64, "((H:3,(C:1,B:1):2):1,G:4);");

julia> adjacency_tree(ptree)
Directed graph with 7 nodes and the following labeling(s):
label: leaves
1 ->
2 ->
3 -> H
4 ->
5 -> C
6 -> B
7 -> G
label: distance
(1, 2) -> 1.0
(1, 7) -> 4.0
(2, 3) -> 3.0
(2, 4) -> 2.0
(4, 5) -> 1.0
(4, 6) -> 1.0
```
"""
function adjacency_tree(ptree::PhylogeneticTree{T}; add_labels=true) where T
  udir_tree = Graph{Undirected}(ptree.pm_ptree.ADJACENCY)
  n = nv(udir_tree)
  dir_tree = Graph{Directed}(n)
  edge_lengths = pm_object(ptree).EDGE_LENGTHS
  edge_labels = Dict{NTuple{2, Int}, T}()

  queue = [1]
  visited = fill(false, n)
  visited[1] = true
  while length(queue) > 0
    x = popfirst!(queue)
    for y in neighbors(udir_tree, x)
      node = vertex_perm(ptree)[y]
      if visited[y] == false
        add_edge!(dir_tree, vertex_perm(ptree)[x], node)
        edge_labels[vertex_perm(ptree)[x], node] = edge_lengths[x, y]
        push!(queue, y)
        visited[y] = true
      end
    end
  end
  add_labels && label!(dir_tree, edge_labels, nothing; name=:distance)
  add_labels && label!(dir_tree,
                       nothing,
                       Dict{Int, String}(
                         v => pm_object(ptree).LABELS[vertex_perm(ptree)[v]] for v in 1:n ); name=:leaves)
  return dir_tree
end

# not exported yet
root(pt::PhylogeneticTree) = root(adjacency_tree(pt))

function Base.show(io::IO, m::MIME"text/plain", ptree::PhylogeneticTree{T}) where T
  print(io, "Phylogenetic tree with $T type coefficients")
end
  
function Base.show(io::IO, ptree::PhylogeneticTree{T}) where T
  print(io, "Phylogenetic tree with $T type coefficients")
end

@doc raw"""
    is_equidistant(ptree::PhylogeneticTree)

Check if the phylogenetic tree `ptree` is equidistant.

# Examples
Make a phylogenetic tree with given Newick format and check if it is equidistant.

```jldoctest
julia> ptree = phylogenetic_tree(Float64, "((H:3,(C:1,B:1):2):1,G:4);");

julia> is_equidistant(ptree)
true
```
"""
function is_equidistant(ptree::PhylogeneticTree)
  return pm_object(ptree).EQUIDISTANT::Bool
end


@doc raw"""
    cophenetic_matrix(ptree::PhylogeneticTree)

Return the cophenetic matrix of the phylogenetic tree `ptree`.

# Examples
Make a phylogenetic tree with given Newick format and print its cophenetic matrix.

```jldoctest
julia> ptree = phylogenetic_tree(Float64, "((H:3,(C:1,B:1):2):1,G:4);");

julia> cophenetic_matrix(ptree)
4×4 Matrix{Float64}:
 0.0  2.0  8.0  6.0
 2.0  0.0  8.0  6.0
 8.0  8.0  0.0  8.0
 6.0  6.0  8.0  0.0
```
"""
function cophenetic_matrix(ptree::PhylogeneticTree{Float64})
  return convert(Matrix, pm_object(ptree).COPHENETIC_MATRIX)::Matrix{Float64}
end

function cophenetic_matrix(ptree::PhylogeneticTree{QQFieldElem})
  return matrix(QQ, pm_object(ptree).COPHENETIC_MATRIX)::QQMatrix
end

@doc raw"""
    taxa(ptree::PhylogeneticTree)

Return the taxa of the phylogenetic tree `ptree`.

# Examples
Make a phylogenetic tree with given Newick format and print its taxa.

```jldoctest
julia> ptree = phylogenetic_tree(Float64, "((H:3,(C:1,B:1):2):1,G:4);");

julia> taxa(ptree)
4-element Vector{String}:
 "B"
 "C"
 "G"
 "H"
```
"""
function taxa(ptree::PhylogeneticTree)
  return convert(Array{String}, pm_object(ptree).TAXA)::Array{String}
end

@doc raw"""
    newick(ptree::PhylogeneticTree)

Return a Newick representation of the phylogenetic tree `ptree`.

# Examples
Make a phylogenetic tree from a matrix and print a Newick representation of it.

```jldoctest
julia> mat = [0. 2 8 6; 2 0 8 6; 8 8 0 8; 6 6 8 0]
4×4 Matrix{Float64}:
 0.0  2.0  8.0  6.0
 2.0  0.0  8.0  6.0
 8.0  8.0  0.0  8.0
 6.0  6.0  8.0  0.0

julia> tax = ["Bonobo", "Chimpanzee", "Gorilla", "Human"]
4-element Vector{String}:
 "Bonobo"
 "Chimpanzee"
 "Gorilla"
 "Human"

julia> tree_mat = phylogenetic_tree(mat, tax);

julia> newick(tree_mat)
"Gorilla:4,(Human:3,(Bonobo:1,Chimpanzee:1):2):1;"
```
"""
function newick(ptree::PhylogeneticTree)
  return convert(String, pm_object(ptree).NEWICK)::String
end


@doc raw"""
    tropical_median_consensus(arr::Vector{PhylogeneticTree{T}})

Compute the tropical median consensus tree of the equidistant
phylogenetic trees from the vector `arr`.  The output is then
equidistant, too.  The method is robust (ie., the consensus tree
dpends continuosuly on the input), and it is Pareto and co-Pareto on
rooted triplets.

The algorithm, based on tropical convexity and optimal transport, is
explained in [CJ24](@cite).

# Examples
Compute the tropical median consensus of three trees and print one of its
Newick representations.

```jldoctest
julia> t1 = phylogenetic_tree(Float64, "((H:30,(C:10,B:10):20):10,G:40);");

julia> t2 = phylogenetic_tree(Float64, "(((H:10,C:10):20,B:30):10,G:40);");

julia> t3 = phylogenetic_tree(Float64, "((H:25,C:25):15,(B:15,G:15):25);");

julia> arr = [t1, t2, t3];

julia> tc = tropical_median_consensus(arr);

julia> newick(tc)
"G:40,(B:35,(C:30,H:30):5):5;"
```
"""
function tropical_median_consensus(arr::Vector{PhylogeneticTree{T}}) where {T <: Union{Float64, QQFieldElem}}

  n = length(arr)
  @req n > 0 "The vector must not be empty"

  phylo_type = Polymake.bigobject_type(pm_object(first(arr)))
  pm_arr = Polymake.Array{Polymake.BigObject}(phylo_type, n)

  pm_arr .= pm_object.(arr)

  pm_cons_tree = Polymake.tropical.tropical_median_consensus(pm_arr)
  return PhylogeneticTree{T}(pm_cons_tree)
end


@doc raw"""
    tropical_median_consensus(trees::Vararg{PhylogeneticTree, N}) where {N}

Compute the tropical median consensus tree of any number of phylogenetic trees
given as parameters.

# Examples
Compute the tropical median consensus of three trees and print one of its
Newick representations.

```jldoctest
julia> t1 = phylogenetic_tree(Float64, "((H:30,(C:10,B:10):20):10,G:40);");

julia> t2 = phylogenetic_tree(Float64, "(((H:10,C:10):20,B:30):10,G:40);");

julia> t3 = phylogenetic_tree(Float64, "((H:25,C:25):15,(B:15,G:15):25);");

julia> tc = tropical_median_consensus(t1, t2, t3);

julia> newick(tc)
"G:40,(B:35,(C:30,H:30):5):5;"
```
"""
function tropical_median_consensus(trees::Vararg{PhylogeneticTree, N}) where {N}
  return tropical_median_consensus(collect(trees))
end

@doc raw"""
    leaves(tree::PhylogeneticTree)

Return the indices of the leaves of the `PhylogeneticTree`.

# Examples
```jldoctest
julia> ptree = phylogenetic_tree(Float64, "((H:3,(C:1,B:1):2):1,G:4);");

julia> leaves(ptree)
4-element Vector{Int64}:
 3
 5
 6
 7
```
"""
leaves(pt::PhylogeneticTree) = findall(iszero, outdegree(adjacency_tree(pt)))

interior_nodes(pt::PhylogeneticTree) = findall(>(1), outdegree(adjacency_tree(pt)))
edges(pt::PhylogeneticTree) = edges(adjacency_tree(pt))
n_edges(pt::PhylogeneticTree) = n_edges(adjacency_tree(pt))
is_isomorphic(pt::PhylogeneticTree, g::Graph) = is_isomorphic(adjacency_tree(pt), g)
is_isomorphic(g::Graph, pt::PhylogeneticTree) = is_isomorphic(pt, g)
is_isomorphic(pt1::PhylogeneticTree, pt2::PhylogeneticTree) = is_isomorphic(adjacency_tree(pt1), adjacency_tree(pt2))

@doc raw"""
    n_leaves(tree::PhylogeneticTree)

Return the indices of the leaves of the `PhylogeneticTree`.

# Examples
```jldoctest
julia> ptree = phylogenetic_tree(Float64, "((H:3,(C:1,B:1):2):1,G:4);");

julia> n_leaves(ptree)
4
```
"""
n_leaves(pt::PhylogeneticTree) = n_leaves(adjacency_tree(pt))

