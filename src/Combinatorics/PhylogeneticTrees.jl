struct PhylogeneticTree{T <: Union{Float64, QQFieldElem}}
  pm_ptree::Polymake.LibPolymake.BigObjectAllocated
end

function pm_object(PT::PhylogeneticTree)
  return PT.pm_ptree
end

@doc raw"""
    phylogenetic_tree(T::Type{<:Union{Float64, QQFieldElem}}, newick::String)

Constructs a phylogenetic tree with Newick representation `newick`. `T` indicates 
the numerical type of the edge lengths.

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
```
"""
function phylogenetic_tree(T::Type{<:Union{Float64, QQFieldElem}}, newick::String)
  pm_ptree = Polymake.graph.PhylogeneticTree{Polymake.convert_to_pm_type(T)}(NEWICK = newick)

  # load graph properties
  pm_ptree.ADJACENCY
  
  return PhylogeneticTree{T}(pm_ptree)
end

@doc raw"""
   phylogenetic_tree(M::Matrix{T}, taxa::Vector{String}) where T <: Union{Float64, QQFieldElem}

Constructs a phylogenetic tree with cophenetic matrix `M` and taxa `taxa`. The matrix `M` must be
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
    adjacency_tree(ptree::PhylogeneticTree)

Returns the underlying graph of the phylogenetic tree `ptree`.

# Examples
Make a phylogenetic tree with given Newick format and print its underlying graph.

```jldoctest
julia> ptree = phylogenetic_tree(Float64, "((H:3,(C:1,B:1):2):1,G:4);");

julia> adjacency_tree(ptree)
Undirected graph with 7 nodes and the following edges:
(2, 1)(3, 2)(4, 2)(5, 4)(6, 4)(7, 1)
```
"""
function adjacency_tree(ptree::PhylogeneticTree)
  return Graph{Undirected}(pm_object(ptree).ADJACENCY)
end

function Base.show(io::IO, ptree::PhylogeneticTree{T}) where T
  print(io, "Phylogenetic tree with $T type coefficients")
end

@doc raw"""
    equidistant(ptree::PhylogeneticTree)

Checks if the phylogenetic tree `ptree` is equidistant.

# Examples
Make a phylogenetic tree with given Newick format and check if it is equidistant.

```jldoctest
julia> ptree = phylogenetic_tree(Float64, "((H:3,(C:1,B:1):2):1,G:4);");

julia> equidistant(ptree)
true
```
"""
function equidistant(ptree::PhylogeneticTree)
  return pm_object(ptree).EQUIDISTANT::Bool
end


@doc raw"""
    cophenetic_matrix(ptree::PhylogeneticTree)

Returns the cophenetic matrix of the phylogenetic tree `ptree`.

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

Returns the taxa of the phylogenetic tree `ptree`.

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

Returns a Newick representation of the phylogenetic tree `ptree`.

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

Computes the tropical median consensus tree of the phylogenetic trees from
the vector `arr`.

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

Computes the tropical median consensus tree of any number of phylogenetic trees
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
