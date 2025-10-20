################################################
#### PHYLOGENETIC DATA STRUCTURES & METHODS ####
################################################


###################################################################################
#
#       Phylogenetic Networks
#
###################################################################################
@doc raw"""
    PhylogeneticNetwork{N, L} <: AbstractGraph{Directed}

A data structure representing a phylogenetic network.

A phylogenetic network is a directed acyclic graph where leaves (nodes with out-degree 0) represent observed species, and internal nodes represent ancestral species. Hybrid nodes (nodes with in-degree greater than 1) represent reticulation events like hybridization.

**Type Parameters:**
* `N`: The number of hybrid nodes in the network.
* `L`: The level of the phylogenetic network.

**Fields:**
* `graph::Graph{Directed}`: The underlying directed graph representing the network topology.
* `hybrids::Dict{Int, Vector{Edge}}`: A dictionary mapping each hybrid vertex ID to its incoming hybrid edges.
"""
struct PhylogeneticNetwork{N, L} <: AbstractGraph{Directed}
  graph::Graph{Directed}
  hybrids::Dict{Int, Vector{Edge}}

  function PhylogeneticNetwork(G::Graph{Directed},
                              hybrid_edgs::Union{Nothing, Vector{Edge}} = nothing)

    # TODO: If hybrid_edgs is not empty, check that they are correct. Maybe the input should just be the graph 
    if isnothing(hybrid_edgs)
      h_nodes = hybrid_vertices(G)
      hybrid_edgs = Dict{Int, Vector{Edge}}(i => hybrid_edges(G, i) for i in h_nodes)
    end

    if !is_phylogenetic_network(G)
      error("The graph `$(collect(edges(G)))` 
        with $(length(hybrid_edgs)) hybrid nodes is not a phylogenetic network.")
    end

    level = level_phylogenetic_network(G)
    return new{length(hybrid_edgs), level}(G, hybrid_edgs)
  end
end

@doc raw"""
    phylogenetic_network(G::Graph{Directed})

A constructor for creating a `PhylogeneticNetwork` from a directed graph `G`.
It automatically detects hybrid nodes and edges.

# Example
```jldoctest
julia> G1 = graph_from_edges(Directed,[[5,6], [5,4], [6,4], [5,2], [6,3], [4,1]]);

julia> phylogenetic_network(G1)
Level-1 phylogenetic network with hybrid nodes {4} and edges
  (4, 1)(5, 2)(5, 4)(5, 6)(6, 3)(6, 4)

julia> G2 = graph_from_edges(Directed,[[7,2], [8,7], [7,6], [8,6], [9,8], [10,9], [11,10], [11,12], [12,10], [7,2], [6,1], [9,3], [11,4], [12,5]]);

julia> phylogenetic_network(G2)
Level-1 phylogenetic network with hybrid nodes {6, 10} and edges
  (6, 1)(7, 2)(7, 6)(8, 6)(8, 7)(9, 3)(9, 8)(10, 9)(11, 4)(11, 10)(11, 12)(12, 5)(12, 10)

julia> G3 = graph_from_edges(Directed,[[6,2], [6,7], [6,5], [7,5], [5,4], [7,8], [8,4], [4,1], [8,3]]);

julia> phylogenetic_network(G3)
Level-2 phylogenetic network with hybrid nodes {4, 5} and edges
  (4, 1)(5, 4)(6, 2)(6, 5)(6, 7)(7, 5)(7, 8)(8, 3)(8, 4)
```
"""
function phylogenetic_network(G::Graph{Directed})
    return PhylogeneticNetwork(G, nothing)
end

function Base.show(io::IO,  m::MIME"text/plain", N::PhylogeneticNetwork)
  l = level(N)
  E = join(["($(src(e)), $(dst(e)))" for e in edges(N)])
  hybN = join(hybrid_vertices(N), ", ")
  print(io, "Level-$l phylogenetic network with hybrid nodes {$hybN} and edges
  $E")  
end

@doc raw"""
    level(::PhylogeneticNetwork{N, L}) where {N, L}

Return the level `L` of the phylogenetic network.

# Example
```jldoctest
julia> N = Oscar.phylogenetic_network(graph_from_edges(Directed,[[5,6], [5,4], [6,4], [5,2], [6,3], [4,1]]));

julia> level(N)
1
```
"""
level(::PhylogeneticNetwork{N, L}) where {N, L} = L

@doc raw"""
    n_hybrid(::PhylogeneticNetwork{N, L}) where {N, L}

Return the number of hybrid nodes `N` in the phylogenetic network.

# Example
```jldoctest
julia> N = Oscar.phylogenetic_network(graph_from_edges(Directed,[[5,6], [5,4], [6,4], [5,2], [6,3], [4,1]]));

julia> n_hybrid(N)
1 
```
"""
n_hybrid(::PhylogeneticNetwork{N, L}) where {N, L} = N

@doc raw"""
    graph(N::PhylogeneticNetwork)

Return the underlying `Graph{Directed}` of the phylogenetic network.

# Example
```jldoctest
julia> N = Oscar.phylogenetic_network(graph_from_edges(Directed,[[5,6], [5,4], [6,4], [5,2], [6,3], [4,1]]));

julia> graph(N)
Directed graph with 6 nodes and the following edges:
(4, 1)(5, 2)(5, 4)(5, 6)(6, 3)(6, 4)
```
"""
graph(N::PhylogeneticNetwork) = N.graph

@doc raw"""
    hybrids(N::PhylogeneticNetwork)

Return the dictionary mapping hybrid vertex IDs to their incoming hybrid edges.

# Example
```jldoctest
julia> N = Oscar.phylogenetic_network(graph_from_edges(Directed,[[5,6], [5,4], [6,4], [5,2], [6,3], [4,1]]));

julia> hybrids(N)
Dict{Int64, Vector{Edge}} with 1 entry:
  4 => [Edge(5, 4), Edge(6, 4)]
```
"""
hybrids(N::PhylogeneticNetwork) = N.hybrids
@doc raw"""
    hybrid_vertices(N::PhylogeneticNetwork)

Return a list of the hybrid vertices in the phylogenetic network.

# Example
```jldoctest
julia> N = Oscar.phylogenetic_network(graph_from_edges(Directed,[[5,6], [5,4], [6,4], [5,2], [6,3], [4,1]]));

julia> hybrid_vertices(N)
1-element Vector{Int64}:
 4
```
"""
hybrid_vertices(N::PhylogeneticNetwork) = hybrid_vertices(graph(N))

@doc raw"""
    hybrid_edges(N::PhylogeneticNetwork)

Return a list of all hybrid edges in the phylogenetic network.

# Example
```jldoctest
julia> N = Oscar.phylogenetic_network(graph_from_edges(Directed,[[5,6], [5,4], [6,4], [5,2], [6,3], [4,1]]));

julia> hybrid_edges(N)
1-element Vector{Vector{Edge}}:
 [Edge(5, 4), Edge(6, 4)]
```
"""
hybrid_edges(N::PhylogeneticNetwork) = hybrid_edges(graph(N))

@doc raw"""
    n_edges(N::PhylogeneticNetwork)

Return the total number of edges in the phylogenetic network.

# Example
```jldoctest
julia> N = Oscar.phylogenetic_network(graph_from_edges(Directed,[[5,6], [5,4], [6,4], [5,2], [6,3], [4,1]]));

julia> n_edges(N)
6
```
"""
n_edges(N::PhylogeneticNetwork) = n_edges(graph(N))

@doc raw"""
    n_leaves(N::PhylogeneticNetwork)

Return the number of leaves in the phylogenetic network.

# Example
```jldoctest
julia> N = Oscar.phylogenetic_network(graph_from_edges(Directed,[[5,6], [5,4], [6,4], [5,2], [6,3], [4,1]]));

julia> n_leaves(N)
3
```
"""
n_leaves(N::PhylogeneticNetwork) = n_leaves(graph(N))

@doc raw"""
    leaves(N::PhylogeneticNetwork)

Return the leaf nodes of the phylogenetic network.

# Example
```jldoctest
julia> N = Oscar.phylogenetic_network(graph_from_edges(Directed,[[5,6], [5,4], [6,4], [5,2], [6,3], [4,1]]));

julia> leaves(N)
3-element Vector{Int64}:
 1
 2
 3
```
"""
leaves(N::PhylogeneticNetwork) = leaves(graph(N))

@doc raw"""
    edges(N::PhylogeneticNetwork)

Return an iterator over the edges of the phylogenetic network.

# Example
```jldoctest
julia> N = Oscar.phylogenetic_network(graph_from_edges(Directed,[[5,6], [5,4], [6,4], [5,2], [6,3], [4,1]]));

julia> collect(edges(N))
6-element Vector{Edge}:
 Edge(4, 1)
 Edge(5, 2)
 Edge(5, 4)
 Edge(5, 6)
 Edge(6, 3)
 Edge(6, 4)
```
"""
edges(N::PhylogeneticNetwork) = edges(graph(N))

###################################################################################
#
#       Structs PhylogeneticModel & GroupBasedPhylogeneticModel
#
###################################################################################

### PhylogeneticModel

@doc raw"""
    PhylogeneticModel{GT, L, M, R}
    
A data structure representing a general phylogenetic model on a directed tree.

This model defines a probability distribution on the states of the leaves of a phylogenetic tree,
based on a root distribution and transition matrices on the edges. The model is defined over a
base field (e.g., `QQ` for rational numbers).

**Type Parameters:**
* `GT`: The graph type, expected to be `Graph{Directed}`, `PhylogeneticTree` or `PhylogeneticNetwork`.
* `L`: The type of the graph labelings, a `NamedTuple`.
* `M`: The type of elements in the transition matrix structure, typically a `VarName` (a `Symbol`) or a polynomial ring element (`MPolyRingElem`).
* `R`: The type of elements in the root distribution vector.

**Fields:**
* `base_field::Field`: The base field of the model's parameters.
* `graph::GT`: The underlying directed graph (tree) of the model.
* `labelings::L`: A named tuple of graph maps for internal use.
* `trans_matrix_structure::Matrix{M}`: A matrix representing the symbolic structure of the transition matrices for each edge.
* `root_distribution::Vector{T}`: A vector representing the distribution at the root of the tree.
* `n_states::Int`: The number of states in the model.
* `model_parameter_name::VarName`: The symbolic name for the variables of the model ring.
"""
@attributes mutable struct PhylogeneticModel{GT, L, M, R} <: GraphicalModel{GT, L}
  base_field::Field
  graph::GT
  labelings::L ## Do we need this?
  trans_matrix_structure::Matrix{M}
  root_distribution::Vector{R}
  n_states::Int
  model_parameter_name::VarName

  # Constructor for PhylogeneticTree
  function PhylogeneticModel(F::Field,
                             G::T, # Is it ok to just have  G <: PhylogeneticTree ?
                             trans_matrix_structure::Matrix,
                             root_distribution::Union{Nothing, Vector} = nothing,
                             varname::VarName="p") where T <: PhylogeneticTree

    n_states = size(trans_matrix_structure)[1]
    if isnothing(root_distribution)
      root_distribution = F.(repeat([1//n_states], outer = n_states))
    end
    graph_maps = NamedTuple(_graph_maps(adjacency_tree(G; add_labels=false)))
    graph_maps = isempty(graph_maps) ? nothing : graph_maps

    return new{
      typeof(G),
      typeof(graph_maps),
      typeof(first(trans_matrix_structure)), 
      typeof(first(root_distribution))
    }(F, G, graph_maps, trans_matrix_structure, root_distribution, n_states, varname)
  end

  # Constructor for PhylogeneticNetwork
  function PhylogeneticModel(F::Field,
                             G::N,
                             trans_matrix_structure::Matrix,
                             root_distribution::Union{Nothing, Vector} = nothing,
                             varnames::VarName="p") where N <: PhylogeneticNetwork

    if n_hybrid(G) == 0
      @warn("The phylogenetic network has no hybrid nodes and is therefore a phylogenetic tree. Consider converting it to a PhylogeneticTree for efficiency.")
    end

    n_states = size(trans_matrix_structure)[1]
    if isnothing(root_distribution)
      root_distribution = F.(repeat([1//n_states], outer = n_states))
    end

    graph_maps = NamedTuple(_graph_maps(graph(G)))
    graph_maps = isempty(graph_maps) ? nothing : graph_maps
    return new{typeof(G), typeof(graph_maps),
               typeof(first(trans_matrix_structure)), 
               typeof(first(root_distribution))}(F, G,
                                                 graph_maps,
                                                 trans_matrix_structure,
                                                 root_distribution,
                                                 n_states,
                                                 varnames)
  end

  # Constructor for Graph{Directed}
  function PhylogeneticModel(F::Field,
                             G::Graph{Directed},
                             trans_matrix_structure::Matrix,
                             root_distribution::Union{Nothing, Vector} = nothing,
                             varname::VarName="p")

    if _is_tree(G)
      pt = phylogenetic_tree(QQFieldElem, G)
      return PhylogeneticModel(F, pt, trans_matrix_structure, root_distribution, varname)
    end

    if is_phylogenetic_network(G)
      N = phylogenetic_network(G)
      return PhylogeneticModel(F, N, trans_matrix_structure, root_distribution, varname)
    end

    error("$G is not a phylogenetic tree nor a phylogenetic network")
  end

  function PhylogeneticModel(G::AbstractGraph{Directed},
                             trans_matrix_structure::Matrix{M},
                             root_distribution::Vector{T},
                             varnames::VarName="p") where {M <:  MPolyRing{<:MPolyRingElem}, T <: MPolyRing}
    trans_ring = parent(first(trans_matrix_structure))
    root_ring = base_ring(trans_ring)
    F = coefficient_ring(root_ring)
    
    if root_ring != parent(first(root_distribution))
      error("Rings for the root distribution parameters do not match. Expected ring `$(root_ring)` but got `$(parent(first(root_distribution)))`.")
    end

    PhylogeneticModel(F, G, trans_matrix_structure, root_distribution, varnames)
  end

  function PhylogeneticModel(G::AbstractGraph{Directed},
                             trans_matrix_structure::Matrix,
                             root_distribution::Union{Nothing, Vector} = nothing,
                             varnames::VarName="p")
    return PhylogeneticModel(QQ, G, trans_matrix_structure, 
                             root_distribution, varnames)
  end

end

@doc raw"""
    phylogenetic_model(F::Field, G::AbstractGraph{Directed}, trans_matrix_structure::Matrix, root_distribution::Union{Nothing, Vector} = nothing, varnames::VarName="p")
    phylogenetic_model(G::AbstractGraph{Directed}, trans_matrix_structure::Matrix, root_distribution::Union{Nothing, Vector} = nothing, varnames::VarName="p")
    phylogenetic_model(G::AbstractGraph{Directed}, trans_matrix_structure::Matrix{M}, root_distribution::Vector{T}, varnames::VarName="p") where {M <: MPolyRing{<:MPolyRingElem}, R <: MPolyRing}
      
Construct a `PhylogeneticModel` from a field `F`, a graph `G`, a symbolic transition matrix `trans_matrix_structure`, and an optional `root_distribution`.

`G` can either be of type `Graph{Directed}`, `PhylogeneticTree` or `PhylogeneticNetwork`.

If `root_distribution` is not provided, a uniform distribution is assumed. If `F` is not provided, it constructs a `PhylogeneticModel` using the default rational field (`QQ`).

If `trans_matrix_structure` is a `Matrix{<: MPolyRing{<:MPolyRingElem}}` and `root_distribution` is a `Vector{<: MPolyRing}`, then the constructor ensures the rings of the root distribution parameters and the transition matrices are compatible.

# Example
```jldoctest
julia> tree = phylogenetic_tree(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> M = [:m11 :m12; :m21 :m22];

julia> root_dist = [:r1, :r2];

julia> PM = phylogenetic_model(tree, M, root_dist)
Phylogenetic model on a tree with 3 leaves and 3 edges
with root distribution [r1, r2] and transition matrices of the form
 [:m11 :m12;
  :m21 :m22].
```
"""
function phylogenetic_model(F::Field, G::AbstractGraph{Directed},
                            trans_matrix_structure::Matrix,
                            root_distribution::Union{Nothing, Vector} = nothing,
                            varnames::VarName="p")
  PhylogeneticModel(F, G, trans_matrix_structure, root_distribution, varnames)
end

function phylogenetic_model(G::AbstractGraph{Directed}, trans_matrix_structure::Matrix{M},
                            root_distribution::Vector{T},
                            varnames::VarName="p") where {M <: MPolyRing{<:MPolyRingElem}, T <: MPolyRing}
  PhylogeneticModel(G, trans_matrix_structure, root_distribution, varnames)
end

function phylogenetic_model(G::AbstractGraph{Directed}, trans_matrix_structure::Matrix,
                            root_distribution::Union{Nothing, Vector} = nothing,
                            varnames::VarName="p")
  return PhylogeneticModel(G, trans_matrix_structure, root_distribution, varnames)
end

### GroupBasedPhylogeneticModel

@doc raw"""
    GroupBasedPhylogeneticModel{GT, L} <: GraphicalModel{GT, L}

A data structure representing a group-based phylogenetic model.

This model is a specialization of `PhylogeneticModel` where the transition matrices and
leaf probabilities can be described by a set of Fourier parameters related to a group structure.

**Type Parameters:**
* `GT`: The graph type, expected to be `Graph{Directed}`.
* `L`: The type of the graph labelings.

**Fields:**
* `phylo_model::PhylogeneticModel`: The underlying standard phylogenetic model.
* `fourier_param_structure::Vector{<: VarName}`: A vector representing the symbolic structure of the Fourier parameters.
* `group::Vector{FinGenAbGroupElem}`: The elements of the finite abelian group used for the model.
* `model_parameter_name::VarName`: The symbolic name for the model's Fourier parameters.
"""
@attributes mutable struct GroupBasedPhylogeneticModel{GT, L} <: GraphicalModel{GT, L}
  phylo_model::PhylogeneticModel{GT, L, <: VarName, <: RingElem} 
  fourier_param_structure::Vector{<: VarName}
  group::Vector{FinGenAbGroupElem}
  model_parameter_name::VarName

  function GroupBasedPhylogeneticModel(pm::PhylogeneticModel{GT, <: VarName, L, <: RingElem},
                                       fourier_param_structure::Vector{<: VarName},
                                       group::Union{Nothing, Vector{FinGenAbGroupElem}} = nothing,
                                       varnames_group_based::VarName="q") where {GT, L}

    if isnothing(group)
      # this feels weird to me, I dont see why you can't just use the group itself?
      # seems like you don't like the ordering?
      group = collect(abelian_group(2,2))
      group = [group[1],group[3],group[2],group[4]]
    end
    new{GT, L}(pm, fourier_param_structure, group, varnames_group_based)
  end
  

  ## I changed this function, so we don't have to define different functions
  ## for phylotree / phylonetwork / directed graph. Is that ok?
  function GroupBasedPhylogeneticModel(F::Field, 
                                       G::AbstractGraph{Directed},
                                       trans_matrix_structure::Matrix{<: VarName},
                                       fourier_param_structure::Vector{<: VarName},
                                       group::Union{Nothing, Vector{FinGenAbGroupElem}} = nothing,
                                       root_distribution::Union{Nothing, Vector} = nothing,
                                       varnames_phylo_model::VarName="p",
                                       varnames_group_based::VarName="q")
    PM = PhylogeneticModel(F, G, trans_matrix_structure, 
                           root_distribution,
                           varnames_phylo_model)
    graph_maps_type = _label_type(PM)
    n = n_states(PM)
    if isnothing(group)
      # this looks odd, can we use a matrix to create the abelian group so we get the order we want
      if n == 4
        group = collect(abelian_group(2,2))
        group = [group[1],group[3],group[2],group[4]]
      elseif n == 2
        group = collect(abelian_group(2))
        group = [group[1],group[2]]
      else
        error("No default group available for a model with $n states. Please provide a vector of group elements.")
      end
    end
    
    return new{typeof(graph(PM)), graph_maps_type}(
      PM, fourier_param_structure, group, varnames_group_based)
  end

  function GroupBasedPhylogeneticModel(G::AbstractGraph{Directed},
                                       trans_matrix_structure::Matrix{<: VarName},
                                       fourier_param_structure::Vector{<: VarName},
                                       group::Union{Nothing, Vector{FinGenAbGroupElem}} = nothing,
                                       root_distribution::Union{Nothing, Vector} = nothing,
                                       varnames_phylo_model::VarName="p",
                                       varnames_group_based::VarName="q")
    return GroupBasedPhylogeneticModel(QQ, G, trans_matrix_structure, fourier_param_structure,
                                       group, root_distribution, 
                                       varnames_phylo_model, varnames_group_based)
  end

end

@doc raw"""
    group_based_phylogenetic_model(F::Field, G::AbstractGraph{Directed}, trans_matrix_structure::Matrix{<: VarName}, fourier_param_structure::Vector{<: VarName}, group::Union{Nothing, Vector{FinGenAbGroupElem}} = nothing, root_distribution::Union{Nothing, Vector} = nothing, varnames_phylo_model::VarName="p", varnames_group_based::VarName="q"))
    group_based_phylogenetic_model(G::AbstractGraph{Directed}, trans_matrix_structure::Matrix{<: VarName}, fourier_param_structure::Vector{<: VarName}, group::Union{Nothing, Vector{FinGenAbGroupElem}} = nothing, root_distribution::Union{Nothing, Vector} = nothing, varnames_phylo_model::VarName="p", varnames_group_based::VarName="q")

Construct a `GroupBasedPhylogeneticModel` from a field `F`, a graph `G` (of type `Graph{Directed}`, `PhylogeneticTree` or `PhylogeneticNetwork`), a symbolic transition matrix, symbolic Fourier parameters, an optional group, and optional root distribution. 

If the `group` is not provided, a default is used based on the number of states: the group $Z_2$ for 2-state models, and the Klein four-group ($Z_2 \times Z_2$) for 4-state models. For other numbers of states, a group must be specified.
If the `root_distribution` is not provided, a uniform distribution is assumed. If `F` is not provided, the model is constructed over the default rational field (`QQ`).

# Example
```jldoctest
julia> tree = phylogenetic_tree(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> M = [:a :b; :b :a];

julia> fourier_params = [:x, :y];

julia> PM = group_based_phylogenetic_model(tree, M, fourier_params)
Group-based phylogenetic model on a tree with 3 leaves and 3 edges
with root distribution [1//2, 1//2],
transition matrices of the form
 [:a :b;
  :b :a]
and fourier parameters of the form [:x, :y].
```
"""
function group_based_phylogenetic_model(F::Field, 
                                      G::AbstractGraph{Directed},
                                      trans_matrix_structure::Matrix{<: VarName},
                                      fourier_param_structure::Vector{<: VarName},
                                      group::Union{Nothing, Vector{FinGenAbGroupElem}} = nothing,
                                      root_distribution::Union{Nothing, Vector} = nothing,
                                      varnames_phylo_model::VarName="p",
                                      varnames_group_based::VarName="q")
   GroupBasedPhylogeneticModel(F, G, trans_matrix_structure, fourier_param_structure, group,
                               root_distribution, varnames_phylo_model, varnames_group_based)
end

function group_based_phylogenetic_model(G::AbstractGraph{Directed},
                                      trans_matrix_structure::Matrix{<: VarName},
                                      fourier_param_structure::Vector{<: VarName},
                                      group::Union{Nothing, Vector{FinGenAbGroupElem}} = nothing,
                                      root_distribution::Union{Nothing, Vector} = nothing,
                                      varnames_phylo_model::VarName="p",
                                      varnames_group_based::VarName="q")
  return GroupBasedPhylogeneticModel(G, trans_matrix_structure, fourier_param_structure,
                                      group, root_distribution, 
                                      varnames_phylo_model, varnames_group_based)
end


### Properties

_label_type(::PhylogeneticModel{GT, L, M, R}) where {GT, L, M, R} = L

@doc raw"""
    phylogenetic_model(PM::GroupBasedPhylogeneticModel)

Return the underlying `PhylogeneticModel`.

# Example
```jldoctest
julia> tree = phylogenetic_tree(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> PM = jukes_cantor_model(tree);

julia> phylogenetic_model(PM)
Phylogenetic model on a tree with 3 leaves and 3 edges
with root distribution [1//4, 1//4, 1//4, 1//4] and transition matrices of the form
 [:a :b :b :b;
  :b :a :b :b;
  :b :b :a :b;
  :b :b :b :a]. 
```
"""
phylogenetic_model(PM::GroupBasedPhylogeneticModel) = PM.phylo_model


@doc raw"""
    n_states(PM::PhylogeneticModel)
    n_states(PM::GroupBasedPhylogeneticModel)

Return the number of states of the underlying phylogenetic model.

# Example
```jldoctest
julia> tree = phylogenetic_tree(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> n_states(general_markov_model(tree))
4

julia> n_states(cavender_farris_neyman_model(tree))
2
```
"""
n_states(PM::PhylogeneticModel) = PM.n_states

n_states(PM::GroupBasedPhylogeneticModel) = n_states(phylogenetic_model(PM))


@doc raw"""
    transition_matrix(PM::PhylogeneticModel)
    transition_matrix(PM::GroupBasedPhylogeneticModel)

Return the structure of the transition matrices of the underlying phylogenetic model.

# Example
```jldoctest
julia> tree = phylogenetic_tree(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> transition_matrix(general_markov_model(tree))
4×4 Matrix{Symbol}:
 :m11  :m12  :m13  :m14
 :m21  :m22  :m23  :m24
 :m31  :m32  :m33  :m34
 :m41  :m42  :m43  :m44

julia> transition_matrix(cavender_farris_neyman_model(tree))
2×2 Matrix{Symbol}:
 :a  :b  
 :b  :a
```
"""
transition_matrix(PM::PhylogeneticModel) = PM.trans_matrix_structure

transition_matrix(PM::GroupBasedPhylogeneticModel) = transition_matrix(phylogenetic_model(PM))

@doc raw"""
    root_distribution(PM::PhylogeneticModel)
    root_distribution(PM::GroupBasedPhylogeneticModel)

Return the root distribution vector of the underlying phylogenetic model.

# Example
```jldoctest
julia> tree = phylogenetic_tree(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> root_distribution(general_markov_model(tree))
4-element Vector{Symbol}:
 :π1
 :π2
 :π3
 :π4

julia> root_distribution(cavender_farris_neyman_model(tree))
2-element Vector{QQFieldElem}:
 1//2
 1//2
```
"""
root_distribution(PM::PhylogeneticModel) = PM.root_distribution

root_distribution(PM::GroupBasedPhylogeneticModel) = root_distribution(phylogenetic_model(PM))


@doc raw"""
    base_field(PM::PhylogeneticModel)
    base_field(PM::GroupBasedPhylogeneticModel)

Return the base field of the underlying phylogenetic model.

# Example
```jldoctest
julia> tree = phylogenetic_tree(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> base_field(general_markov_model(tree))
Rational field

julia> base_field(jukes_cantor_model(tree))
Rational field
```
"""
base_field(PM::PhylogeneticModel) = PM.base_field

base_field(PM::GroupBasedPhylogeneticModel) = base_field(phylogenetic_model(PM))


@doc raw"""
    varnames(PM::PhylogeneticModel)
    varnames(PM::GroupBasedPhylogeneticModel)

Return the symbolic name for the probability coordinates if `PM` is a `PhylogeneticModel` and for the Fourier corrdinates if `PM` is a `GroupBasedPhylogeneticModel`. 

# Example
```jldoctest
julia> tree = phylogenetic_tree(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> PM = general_markov_model(tree);

julia> varnames(PM)
"p"

julia> varnames(jukes_cantor_model(tree))
"q"
```
"""
varnames(PM::PhylogeneticModel) = PM.model_parameter_name

varnames(PM::GroupBasedPhylogeneticModel) = PM.model_parameter_name

## Do we want to export/document these 2 functions? Do we even need them?
varnames_probabilities(PM::GroupBasedPhylogeneticModel) = varnames(phylogenetic_model(PM))

varnames_fourier(PM::GroupBasedPhylogeneticModel) = varnames(PM)

@doc raw"""
    graph(PM::GroupBasedPhylogeneticModel)

Return the graph of the underlying phylogenetic model.

# Example
```jldoctest
julia> tree = phylogenetic_tree(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> PM = jukes_cantor_model(tree);

julia> graph(PM)
Phylogenetic tree with QQFieldElem type coefficients
```
"""
graph(PM::GroupBasedPhylogeneticModel) = graph(phylogenetic_model(PM))

@doc raw"""
    group(PM::GroupBasedPhylogeneticModel)

Return the finite abelian group associated with the model.

# Example
```jldoctest
julia> tree = phylogenetic_tree(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> PM = jukes_cantor_model(tree);

julia> group(PM)
4-element Vector{FinGenAbGroupElem}:
 [0, 0]
 [0, 1]
 [1, 0]
 [1, 1]
```
"""
group(PM::GroupBasedPhylogeneticModel) = PM.group

@doc raw"""
    fourier_parameters(PM::GroupBasedPhylogeneticModel)

Return the symbolic structure of the Fourier parameters.

# Example
```jldoctest
julia> tree = phylogenetic_tree(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> PM = jukes_cantor_model(tree);

julia> fourier_parameters(PM)
4-element Vector{Symbol}:
 :x
 :y
 :y
 :y
```
"""
fourier_parameters(PM::GroupBasedPhylogeneticModel) = PM.fourier_param_structure

function Base.show(io::IO, PM::PhylogeneticModel{<:PhylogeneticTree})
  gr = graph(PM)

  nl = length(leaves(gr))
  ne = length(collect(edges(gr)))
  root_dist = join(root_distribution(PM), ", " )
 
  print(io, "Phylogenetic model on a tree with $(nl) leaves and $(ne) edges \n") # \n )
  print(io, "with root distribution [$(root_dist)] ")
  print(io, "and transition matrices of the form \n ")

  M = string(PM.trans_matrix_structure)
  M = split(M, "[", limit=2)[2]
  print(io, "[", replace(M, ";" => ";\n "))
  print(io, ". ")
end

function Base.show(io::IO, PM::PhylogeneticModel{<:PhylogeneticNetwork})
  gr = graph(PM)

  nl = length(leaves(gr))
  ne = length(collect(edges(gr)))
  root_dist = join(root_distribution(PM), ", " )

  print(io, "Phylogenetic model on a level-$(level(gr)) network with $(n_hybrid(gr)) hybrid node, ")
  print(io, "$(nl) leaves and $(ne) edges \n") # \n )
  print(io, "with root distribution [$(root_dist)] ")
  print(io, "and transition matrices of the form \n ")

  M = string(PM.trans_matrix_structure)
  M = split(M, "[", limit=2)[2]
  print(io, "[", replace(M, ";" => ";\n "))
  print(io, ". ")
end

function Base.show(io::IO, PM::GroupBasedPhylogeneticModel{<: PhylogeneticTree})

  gr = graph(PM)

  nl = length(leaves(gr))
  ne = length(collect(edges(gr)))
  root_dist = join(PM.phylo_model.root_distribution, ", " )
 
  print(io, "Group-based phylogenetic model on a tree with $(nl) leaves and $(ne) edges \n") # \n )
  print(io, "with root distribution [$(root_dist)], \n")
  print(io, "transition matrices of the form \n ")

  M = string(PM.phylo_model.trans_matrix_structure)
  print(io, replace(M, ";" => ";\n "))
  print(io, "\n")

  print(io, "and fourier parameters of the form ")
  print(io, "$(PM.fourier_param_structure).")

end

function Base.show(io::IO, PM::GroupBasedPhylogeneticModel{<: PhylogeneticNetwork})

  gr = graph(PM)

  nl = length(leaves(gr))
  ne = length(collect(edges(gr)))
  root_dist = join(PM.phylo_model.root_distribution, ", " )
 
  print(io, "Group-based phylogenetic model on a level-$(level(gr)) network with $(n_hybrid(gr)) hybrid node, ")
  print(io, "$(nl) leaves  \nand $(ne) edges ") # \n )
  print(io, "with root distribution [$(root_dist)], \n")
  print(io, "transition matrices of the form \n ")

  M = string(PM.phylo_model.trans_matrix_structure)
  print(io, replace(M, ";" => ";\n "))
  print(io, "\n")

  print(io, "and fourier parameters of the form ")
  print(io, "$(PM.fourier_param_structure).")

end

###################################################################################
#
#       Parameter Rings
#
###################################################################################

### PhylogeneticModel & PhylogeneticTree

#TODO: add explanation sorted_edges in the docs
@doc raw"""
    parameter_ring(PM::PhylogeneticModel; cached=false, sorted_edges::Union{Vector{Edge}, Nothing} = nothing)
   
Create the polynomial ring for the parameters of a phylogenetic model `PhylogeneticModel{GT, L, M, T}` where
`GT <: PhylogeneticTree` or `GT <: PhylogeneticNetwork`.

If `M  <: VarName`, returns a tuple containing:
  1.  The polynomial ring.
  2.  A dictionary mapping parameter variables and edges to the corresponding ring generators.
  3.  The root distribution vector.

If `M  <: MPolyRingElem`, then returns a tuple containing:
  1.  The polynomial ring.
  2.  A dictionary mapping each edge to a homomorphism (`Oscar.MPolyAnyMap`) from the original
      transition matrix ring to the new parameter ring.
  3.  The root distribution vector.

If `GT <: PhylogeneticNetwork` it additionally returns:
  4. A dictionary mapping hybrid edges to the corresponding hybrid parameter in the ring.
"""
@attr Tuple{
  MPolyRing, 
  GraphTransDict,
  Vector{RT}
} function parameter_ring(PM::PhylogeneticModel{<:PhylogeneticTree, L, <: VarName, RT}; cached=false, 
                          sorted_edges::Union{Vector{Edge}, Nothing} = nothing) where {L, RT <: FieldElem} 
  vars = unique(transition_matrix(PM))
  edge_gens = [x => 1:n_edges(graph(PM)) for x in vars]
  R, x... = polynomial_ring(base_field(PM), edge_gens...; cached=cached)
  
  R, Dict{Tuple{VarName, Edge}, MPolyRingElem}(
    (vars[i], e) => x[i][j] for i in 1:length(vars), 
      (j,e) in enumerate(sort_edges(graph(PM), sorted_edges))
      ), root_distribution(PM)
end

@attr Tuple{
  MPolyRing,
  GraphTransDict,
  Vector{<:MPolyRingElem}
} function parameter_ring(PM::PhylogeneticModel{<:PhylogeneticTree, L, <: VarName}; cached=false,
                          sorted_edges::Union{Vector{Edge}, Nothing} = nothing) where {L}
  vars = unique(transition_matrix(PM))
  edge_gens = [x => 1:n_edges(graph(PM)) for x in vars]
  R, r, x... = polynomial_ring(base_field(PM),
                               root_distribution(PM),
                               edge_gens...; cached=cached)

  R, Dict{Tuple{VarName, Edge}, MPolyRingElem}(
    (vars[i], e) => x[i][j] for i in 1:length(vars), 
      (j,e) in enumerate(sort_edges(graph(PM), sorted_edges))
      ), r
end

@attr Tuple{
  MPolyRing, 
  Dict{Edge, Oscar.MPolyAnyMap}, 
  Vector{RT}
} function parameter_ring(PM::PhylogeneticModel{<:PhylogeneticTree, L, <: MPolyRingElem, RT}; cached=false,
                          sorted_edges::Union{Vector{Edge}, Nothing} = nothing) where {L, RT <: FieldElem}
  trans_ring = parent(first(transition_matrix(PM)))
  transition_vars = gens(trans_ring)

  edge_gens = [x => 1:n_edges(graph(PM)) for x in Symbol.(transition_vars)]
  R, x... = polynomial_ring(base_field(PM), edge_gens...; cached=cached)
  
  dict_maps = Dict{Edge, Oscar.MPolyAnyMap}()
  for (j,e) in enumerate(Oscar.sort_edges(graph(PM), sorted_edges))
      map = [x[i][j] for i in 1:length(transition_vars)]
      dict_maps[e] = hom(trans_ring, R, map)
  end

  R, dict_maps, root_distribution(PM)
end

@attr Tuple{
  MPolyRing, 
  Dict{Edge, Oscar.MPolyAnyMap}, 
  Vector{RT}
} function parameter_ring(PM::PhylogeneticModel{<:PhylogeneticTree, L, <: MPolyRingElem, RT}; cached=false,
                          sorted_edges::Union{Vector{Edge}, Nothing} = nothing) where {L, RT <: MPolyRingElem}
  trans_ring = parent(first(transition_matrix(PM)))
  transition_vars = gens(trans_ring)
  root_vars = gens(coefficient_ring(trans_ring))

  edge_gens = [x => 1:n_edges(graph(PM)) for x in Symbol.(transition_vars)]
  R, rv, x... = polynomial_ring(base_field(PM), Symbol.(root_vars), edge_gens..., ; cached=cached)

  coef_map = hom(coefficient_ring(trans_ring), R, rv)

  dict_maps = Dict{Edge, Oscar.MPolyAnyMap}()
  for (j,e) in enumerate(Oscar.sort_edges(graph(PM), sorted_edges))
      map = [x[i][j] for i in 1:length(transition_vars)]
      dict_maps[e] = hom(trans_ring, R, coef_map, map)
  end

  R, dict_maps, rv
end

@attr Tuple{
  MPolyRing{RT}, 
  Dict{Edge, Oscar.MPolyAnyMap}, 
  Vector{RT}
} function parameter_ring(PM::PhylogeneticModel{<:PhylogeneticTree, L, <: MPolyRingElem, RT}; cached=false,
                          sorted_edges::Union{Vector{Edge}, Nothing} = nothing) where {L, RT <: AbstractAlgebra.Generic.RationalFunctionFieldElem} 
  trans_ring = parent(first(transition_matrix(PM)))
  transition_vars = gens(trans_ring)
  root_vars = gens(coefficient_ring(trans_ring))

  edge_gens = [x => 1:n_edges(graph(PM)) for x in Symbol.(transition_vars)]
  R, x... = polynomial_ring(coefficient_ring(trans_ring), edge_gens..., ; cached=cached)

  dict_maps = Dict{Edge, Oscar.MPolyAnyMap}()
  for (j,e) in enumerate(Oscar.sort_edges(graph(PM), sorted_edges))
    map = [x[i][j] for i in 1:length(transition_vars)]
    dict_maps[e] = hom(trans_ring, R, map)
  end
  R, dict_maps, root_vars
end

### PhylogeneticModel & PhylogeneticNetwork
@attr Tuple{
  MPolyRing, 
  GraphTransDict,
  Vector{RT},
  GenDict
} function parameter_ring(PM::PhylogeneticModel{<:PhylogeneticNetwork, L, <: VarName, RT}; cached=false,
                          sorted_edges::Union{Vector{Edge}, Nothing} = nothing) where {L, RT <: FieldElem} 
  N = graph(PM)

  vars = unique(transition_matrix(PM))
  edge_gens = [x => 1:n_edges(N) for x in vars]
  h_nodes = hybrid_vertices(N)

  R, l, x... = polynomial_ring(base_field(PM), :l => (1:length(h_nodes),1:2), edge_gens...; cached=cached)
  
  hyb = hybrids(N)
  R, Dict{Tuple{VarName, Edge}, MPolyRingElem}(
    (vars[i], e) => x[i][j] for i in 1:length(vars), 
      (j,e) in enumerate(sort_edges(N, sorted_edges))
      ), root_distribution(PM), 
      Dict{Edge, MPolyRingElem}(hyb[h_nodes[i]][j] => l[i,j] for i in 1:length(h_nodes) for j in 1:2)
end

@attr Tuple{
  MPolyRing,
  GraphTransDict,
  Vector{<:MPolyRingElem},
  GenDict
} function parameter_ring(PM::PhylogeneticModel{<:PhylogeneticNetwork, L, <: VarName}; cached=false,
                          sorted_edges::Union{Vector{Edge}, Nothing} = nothing)  where {L}
  N = graph(PM)

  vars = unique(transition_matrix(PM))
  edge_gens = [x => 1:n_edges(N) for x in vars]
  h_nodes = hybrid_vertices(N)
  R, r, l, x... = polynomial_ring(base_field(PM),
                               root_distribution(PM),
                               :l => (1:length(h_nodes), 1:2),
                               edge_gens...; cached=cached)

  hyb = hybrids(N)                               
  R, Dict{Tuple{VarName, Edge}, MPolyRingElem}(
    (vars[i], e) => x[i][j] for i in 1:length(vars), 
      (j,e) in enumerate(sort_edges(N, sorted_edges))
      ), r, Dict{Edge, MPolyRingElem}(hyb[h_nodes[i]][j] => l[i,j] for i in 1:length(h_nodes) for j in 1:2)
end

@attr Tuple{
  MPolyRing, 
  Dict{Edge, Oscar.MPolyAnyMap}, 
  Vector{RT},
  GenDict
} function parameter_ring(PM::PhylogeneticModel{<:PhylogeneticNetwork, L, <: MPolyRingElem, RT}; cached=false,
                          sorted_edges::Union{Vector{Edge}, Nothing} = nothing) where {L, RT <: FieldElem}
  N = graph(PM)

  trans_ring = parent(first(transition_matrix(PM)))
  transition_vars = gens(trans_ring)

  edge_gens = [x => 1:n_edges(N) for x in Symbol.(transition_vars)]
  h_nodes = hybrid_vertices(N)
  R, l, x... = polynomial_ring(base_field(PM), :l => (1:length(h_nodes), 1:2), edge_gens...; cached=cached)
  
  dict_maps = Dict{Edge, Oscar.MPolyAnyMap}()
  for (j,e) in enumerate(Oscar.sort_edges(N, sorted_edges))
      map = [x[i][j] for i in 1:length(transition_vars)]
      dict_maps[e] = hom(trans_ring, R, map)
  end

  hyb = hybrids(N)
  R, dict_maps, root_distribution(PM), Dict{Edge, MPolyRingElem}(hyb[h_nodes[i]][j] => l[i,j] for i in 1:length(h_nodes) for j in 1:2)
end

@attr Tuple{
  MPolyRing, 
  Dict{Edge, Oscar.MPolyAnyMap}, 
  Vector{RT},
  GenDict
} function parameter_ring(PM::PhylogeneticModel{<:PhylogeneticNetwork, L, <: MPolyRingElem, RT}; cached=false,
                          sorted_edges::Union{Vector{Edge}, Nothing} = nothing) where {L, RT <: MPolyRingElem}
  N = graph(PM)

  trans_ring = parent(first(transition_matrix(PM)))
  transition_vars = gens(trans_ring)
  root_vars = gens(coefficient_ring(trans_ring))

  edge_gens = [x => 1:n_edges(N) for x in Symbol.(transition_vars)]
  h_nodes = hybrid_vertices(N)
  R, rv, l, x... = polynomial_ring(base_field(PM), Symbol.(root_vars), :l => (1:length(h_nodes), 1:2), edge_gens..., ; cached=cached)

  coef_map = hom(coefficient_ring(trans_ring), R, rv)

  dict_maps = Dict{Edge, Oscar.MPolyAnyMap}()
  for (j,e) in enumerate(Oscar.sort_edges(N, sorted_edges))
      map = [x[i][j] for i in 1:length(transition_vars)]
      dict_maps[e] = hom(trans_ring, R, coef_map, map)
  end

  hyb = hybrids(N)
  R, dict_maps, rv, Dict{Edge, MPolyRingElem}(hyb[h_nodes[i]][j] => l[i,j] for i in 1:length(h_nodes) for j in 1:2)

end

@attr Tuple{
  MPolyRing{R}, 
  Dict{Edge, Oscar.MPolyAnyMap}, 
  Vector{RT},
  GenDict
} function parameter_ring(PM::PhylogeneticModel{<:PhylogeneticNetwork, L, <: MPolyRingElem, RT}; cached=false,
                          sorted_edges::Union{Vector{Edge}, Nothing} = nothing) where {L, RT <: AbstractAlgebra.Generic.RationalFunctionFieldElem}
  N = graph(PM)

  trans_ring = parent(first(transition_matrix(PM)))
  transition_vars = gens(trans_ring)
  root_vars = gens(coefficient_ring(trans_ring))

  edge_gens = [x => 1:n_edges(N) for x in Symbol.(transition_vars)]
  h_nodes = hybrid_vertices(N)
  R, l, x... = polynomial_ring(coefficient_ring(trans_ring), :l => (1:length(h_nodes), 1:2), edge_gens..., ; cached=cached)

  dict_maps = Dict{Edge, Oscar.MPolyAnyMap}()
  for (j,e) in enumerate(Oscar.sort_edges(N, sorted_edges))
    map = [x[i][j] for i in 1:length(transition_vars)]
    dict_maps[e] = hom(trans_ring, R, map)
  end

  hyb = hybrids(N)
  R, dict_maps, root_vars, Dict{Edge, MPolyRingElem}(hyb[h_nodes[i]][j] => l[i,j] for i in 1:length(h_nodes) for j in 1:2)
end

### GroupBasedPhylogeneticModel & PhylogeneticTree
@doc raw"""
    parameter_ring(PM::GroupBasedPhylogeneticModel; cached=false, sorted_edges::Union{Vector{Edge}, Nothing} = nothing)

Create the polynomial ring for the Fourier parameters of the model.

Returns a tuple containing:
1.  The polynomial ring.
2.  A dictionary mapping parameter variables and edges to the corresponding ring generators.

If `GT <: PhylogeneticNetwork` it additionally returns:
  4. A dictionary mapping hybrid edges to the corresponding hybrid parameter in the ring.
"""
@attr Tuple{
  MPolyRing,
  GraphTransDict
} function parameter_ring(PM::GroupBasedPhylogeneticModel{<:PhylogeneticTree}; cached=false,
                          sorted_edges::Union{Vector{Edge}, Nothing} = nothing)
  vars = unique(fourier_parameters(PM))
  edge_gens = [x => 1:n_edges(graph(PM)) for x in vars]
  R, x... = polynomial_ring(base_field(PM), edge_gens...; cached=cached)

  R, Dict{Tuple{VarName, Edge}, MPolyRingElem}(
    (vars[i], e) => x[i][j] for i in 1:length(vars), (j,e) in enumerate(sort_edges(graph(PM), sorted_edges))
  )
end

### GroupBasedPhylogeneticModel & PhylogeneticNetwork
@attr Tuple{
  MPolyRing,
  GraphTransDict,
  GenDict
} function parameter_ring(PM::GroupBasedPhylogeneticModel{<:PhylogeneticNetwork}; cached=false,
                          sorted_edges::Union{Vector{Edge}, Nothing} = nothing)
  N = graph(PM)

  vars = unique(fourier_parameters(PM))
  edge_gens = [x => 1:n_edges(N) for x in vars]
  h_nodes = hybrid_vertices(N)

  R, l, x... = polynomial_ring(base_field(PM), :l => (1:length(h_nodes),1:2), edge_gens...; cached=cached)
  

  hyb = hybrids(N)
  R, Dict{Tuple{VarName, Edge}, MPolyRingElem}(
     (vars[i], e) => x[i][j] for i in 1:length(vars), (j,e) in enumerate(sort_edges(N, sorted_edges))),
     Dict{Edge, MPolyRingElem}(hyb[h_nodes[i]][j] => l[i,j] for i in 1:length(h_nodes) for j in 1:2)

end


###################################################################################
#
#       Model Ring & Parametrizations
#
###################################################################################

@doc raw"""
    full_model_ring(PM::Union{PhylogeneticModel, GroupBasedPhylogeneticModel}; cached=false)

Creates the full model ring in _probability coordinates_ for a `PhylogeneticModel`. Equivalently, creates the full model ring in the _Fourier coordinates_ for a `GroupBasedPhylogeneticModel`.

The output consist on a tuple containing a `ModelRing` and a `Dict` mapping a tuple of states at the leaves to the corresponding coordinate in the ring. This ring has a generator for every possible state configuration at the leaves.

# Example
```jldoctest
julia> tree = phylogenetic_tree(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> R, p = full_model_ring(general_markov_model(tree));

julia> p[1,2,4]
p[1,2,4]
```
"""
@attr Tuple{
  ModelRing{T, U}, 
  Dict{T, U}
} where {T, U} function full_model_ring(PM::Union{PhylogeneticModel, GroupBasedPhylogeneticModel}; cached=false)
  l_indices = leaves_indices(PM)

  return model_ring(base_field(PM), varnames(PM) => l_indices; cached=cached)
end

@doc raw"""
    full_parametrization(PM::PhylogeneticModel{<:PhylogeneticTree})
    full_parametrization(PM::PhylogeneticModel{<:PhylogeneticNetwork})

Constructs the parametrization map from the full model ring in _probability_ coordinates (with generators for all leaf probability configurations) to the parameter ring (with generators from the transition matrices and the root distribution).
"""
@attr MPolyAnyMap function full_parametrization(PM::PhylogeneticModel{<:PhylogeneticTree})
  gr = graph(PM)

  R, _ = full_model_ring(PM)
  S, _ = parameter_ring(PM)

  lvs = sort(leaves(gr))
  lvs_indices = leaves_indices(PM)

  map = [leaves_probability(PM, Dict(lvs[i] => k[i] for i in 1:n_leaves(gr))) for k in lvs_indices]
  hom(R, S, reduce(vcat, map))
end

@attr MPolyAnyMap function full_parametrization(PM::PhylogeneticModel{<:PhylogeneticNetwork})
  N = graph(PM)

  R, _ = full_model_ring(PM)
  S, _ = parameter_ring(PM)

  lvs = sort(leaves(N))
  lvs_indices = leaves_indices(PM)
  h_indices = hybrid_indices(PM)

  t_edges = tree_edges(N)
  hyb = hybrids(N)
  h_nodes = collect(keys(hyb))

  map = [0 for i in lvs_indices]

  for idx in h_indices
    subtree_h_edges = [hyb[h_nodes[i]][idx[i]] for i in eachindex(h_nodes)]
    subtree = graph_from_edges(Directed, vcat(subtree_h_edges, t_edges))

    l = prod([Oscar.entry_hybrid_parameter(PM, e) for e in subtree_h_edges])
    map_subtree = [Oscar.leaves_probability(PM, Dict(lvs[i] => k[i] for i in 1:n_leaves(N)), subtree) for k in lvs_indices]
    map = map + l.*map_subtree
  end

  hom(R, S, reduce(vcat, map))
end

@doc raw"""
    full_parametrization(PM::GroupBasedPhylogeneticModel{<:PhylogeneticTree})
    full_parametrization(PM::GroupBasedPhylogeneticModel{<:PhylogeneticNetwork})

Constructs the parametrization map from the full model ring in _Fourier_ coordinates (with generators for all leaf probability
configurations) to the parameter ring (with generators from the Fourier parameters).
"""
@attr MPolyAnyMap function full_parametrization(PM::GroupBasedPhylogeneticModel{<:PhylogeneticTree})
  gr = graph(PM)

  R, _ = full_model_ring(PM)
  S, _ = parameter_ring(PM)

  lvs = sort(leaves(gr))
  lvs_indices = leaves_indices(PM)

  map = [leaves_fourier(PM, Dict(lvs[i] => k[i] for i in 1:n_leaves(gr))) for k in lvs_indices]
  hom(R, S, reduce(vcat, map))

end

@attr MPolyAnyMap function full_parametrization(PM::GroupBasedPhylogeneticModel{<:PhylogeneticNetwork})
  N = graph(PM)

  R, _ = full_model_ring(PM)
  S, _ = parameter_ring(PM)

  lvs = sort(leaves(N))
  lvs_indices = leaves_indices(PM)
  h_indices = hybrid_indices(PM)

  t_edges = tree_edges(N)
  hyb = hybrids(N)
  h_nodes = collect(keys(hyb))

  map = [0 for i in lvs_indices]

  for idx in h_indices
      subtree_h_edges = [hyb[h_nodes[i]][idx[i]] for i in eachindex(h_nodes)]
      subtree = graph_from_edges(Directed, vcat(subtree_h_edges, t_edges))

      l = prod([Oscar.entry_hybrid_parameter(PM, e) for e in subtree_h_edges])
      map_subtree = [Oscar.leaves_fourier(PM, Dict(lvs[i] => k[i] for i in 1:n_leaves(N)), subtree) for k in lvs_indices]
      map = map + l.*map_subtree
  end

  hom(R, S, reduce(vcat, map))

end

@doc raw"""
    equivalent_classes(PM::Union{PhylogeneticModel, GroupBasedPhylogeneticModel})

Calculates the equivalence classes of the model's parametrizations. Equivalent polynomials
are those that are identical under the chosen parametrization. This is a crucial step for
dimension reduction in the model's algebraic variety.

# Example
```jldoctest
julia> tree = phylogenetic_tree(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> PM = jukes_cantor_model(tree);

julia> R, q = full_model_ring(PM);

julia> class = equivalent_classes(PM);

julia> class[1,2,2]
3-element Vector{MPolyRingElem}:
 q[1,2,2]
 q[1,3,3]
 q[1,4,4]
```
"""
@attr Dict{
      Tuple{Vararg{Int64}}, 
      Vector{MPolyRingElem}
} function equivalent_classes(PM::Union{PhylogeneticModel, GroupBasedPhylogeneticModel}) 
  _, ps = full_model_ring(PM)
  f = full_parametrization(PM);
  
  polys = unique(f.img_gens)
  polys = polys[findall(!is_zero, polys)]

  equivalent_classes = Dict{Tuple{Vararg{Int64}}, Vector{MPolyRingElem}}()
  for poly in polys
      eqv_class = sort([i for i in keys(ps) if f(ps[i]) == poly], rev = false)
      equivalent_classes[eqv_class[1]] = [ps[i...] for i in eqv_class]
  end
  return equivalent_classes
end

@doc raw"""
    model_ring(PM::Union{PhylogeneticModel, GroupBasedPhylogeneticModel}; cached=false)

Creates a reduced model ring based on the equivalence classes of the parametrizations.
This ring has a generator for each unique polynomial (i.e., each equivalence class).

# Example
```jldoctest
julia> tree = phylogenetic_tree(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> R, q = model_ring(cavender_farris_neyman_model(tree));

julia> q[1,1,1]
q[1,1,1]
```

"""
@attr Tuple{
  ModelRing{T, U}, 
  Dict{T, U}
} where {T, U} function model_ring(PM::Union{PhylogeneticModel, GroupBasedPhylogeneticModel}; cached=false)
  eq_classes = equivalent_classes(PM)
  ec_indices = sort(collect(keys(eq_classes)), rev = true)

  return model_ring(base_field(PM), varnames(PM) => ec_indices; cached=cached)
end

@doc raw"""
    parametrization(PM::Union{PhylogeneticModel, GroupBasedPhylogeneticModel})

Constructs the parametrization map from the model ring to the parameter ring.

The rings used depend on the type of `PM`:
- **`PhylogeneticModel`**: Map from probability coordinates to parameters of the transition matrices and root distribution.
- **`GroupBasedPhylogeneticModel`**: Map from Fourier coordinates to Fourier parameters.
"""
function parametrization(PM::Union{PhylogeneticModel, GroupBasedPhylogeneticModel})
  _, p = full_model_ring(PM)
  f = full_parametrization(PM)

  R, _ = model_ring(PM)
  S, _ = parameter_ring(PM)

  lvs_indices = [R.gen_to_index[i] for i in gens(_ring(R))]
  map = [f(p[k...]) for k in lvs_indices]
  hom(R, S, reduce(vcat, map))
end

@doc raw"""
    affine_parametrization(PM::PhylogeneticModel)
    full_affine_parametrization(PM::PhylogeneticModel)

Constructs an affine parametrization (reduced or full) of the model by imposing the transition matrices rows sum to one and the root distribution sums to one.

# Example
```jldoctest
julia> tree = phylogenetic_tree(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> PM = phylogenetic_model(jukes_cantor_model(tree));

julia> affine_parametrization(PM)
Ring homomorphism
  from multivariate polynomial ring in 5 variables over QQ
  to multivariate polynomial ring in 6 variables over QQ
defined by
  p[1,2,3] -> -2*b[1]*b[2]*b[3] + 1//4*b[1]*b[2] + 1//4*b[1]*b[3] + 1//4*b[2]*b[3]
  p[1,2,2] -> 2*b[1]*b[2]*b[3] - 3//4*b[1]*b[2] - 3//4*b[1]*b[3] + 1//4*b[1] + 1//4*b[2]*b[3]
  p[1,2,1] -> 2*b[1]*b[2]*b[3] - 3//4*b[1]*b[2] + 1//4*b[1]*b[3] - 3//4*b[2]*b[3] + 1//4*b[2]
  p[1,1,2] -> 2*b[1]*b[2]*b[3] + 1//4*b[1]*b[2] - 3//4*b[1]*b[3] - 3//4*b[2]*b[3] + 1//4*b[3]
  p[1,1,1] -> -6*b[1]*b[2]*b[3] + 9//4*b[1]*b[2] + 9//4*b[1]*b[3] - 3//4*b[1] + 9//4*b[2]*b[3] - 3//4*b[2] - 3//4*b[3] + 1//4
```
"""
function affine_parametrization(PM::PhylogeneticModel)
  R, p = model_ring(PM) 
  S, _ = parameter_ring(PM)
  f = parametrization(PM)
  ind = [R.gen_to_index[i] for i in gens(_ring(R))]

  edgs = collect(edges(graph(PM)))
  vars = unique([entry_transition_matrix(PM, i, i, e) for i in 1:n_states(PM) for e in edgs])
  vals = unique([(1 - (sum([entry_transition_matrix(PM, i, j, e) for j in 1:n_states(PM)]) - entry_transition_matrix(PM, i, i, e)))  
          for i in 1:n_states(PM) for e in edgs])

  map = [evaluate(f(p[k...]), vars, vals) for k in ind]
  hom(R, S, reduce(vcat, map))
  #TODO: make π sum to one as well. 
end

function full_affine_parametrization(PM::PhylogeneticModel)
  R, p = full_model_ring(PM) 
  S, _ = parameter_ring(PM)
  f = full_parametrization(PM)
  ind = [R.gen_to_index[i] for i in gens(_ring(R))]

  edgs = collect(edges(graph(PM)))
  vars = unique([entry_transition_matrix(PM, i, i, e) for i in 1:n_states(PM) for e in edgs])
  vals = unique([(1 - (sum([entry_transition_matrix(PM, i, j, e) for j in 1:n_states(PM)]) - entry_transition_matrix(PM, i, i, e)))  
          for i in 1:n_states(PM) for e in edgs])

  map = [evaluate(f(p[k...]), vars, vals) for k in ind]
  hom(R, S, reduce(vcat, map))

end

@doc raw"""
    affine_parametrization(PM::GroupBasedPhylogeneticModel)

Constructs an affine parametrization of the group-based model by setting the first
Fourier parameters to 1.

# Example
```jldoctest
julia> tree = phylogenetic_tree(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> PM = jukes_cantor_model(tree);

julia> affine_parametrization(PM)
Ring homomorphism
  from multivariate polynomial ring in 5 variables over QQ
  to multivariate polynomial ring in 6 variables over QQ
defined by
  q[2,3,4] -> y[1]*y[2]*y[3]
  q[2,2,1] -> y[1]*y[2]
  q[2,1,2] -> y[1]*y[3]
  q[1,2,2] -> y[2]*y[3]
  q[1,1,1] -> 1
```
"""
function affine_parametrization(PM::GroupBasedPhylogeneticModel)
  R, q = model_ring(PM)
  S, _ = parameter_ring(PM)
  f = parametrization(PM)
  ind = [R.gen_to_index[i] for i in gens(_ring(R))]

  edgs = collect(edges(graph(PM)))
  vars = [entry_fourier_parameter(PM, 1, e) for e in edgs]
  vals = repeat([1], length(edgs))

  map = [evaluate(f(q[k...]), vars, vals) for k in ind]
  hom(R, S, reduce(vcat, map))
end


###################################################################################
#
#       Fourier - probabilities coordinate change
#
###################################################################################

@doc raw"""
    fourier_transform(PM::GroupBasedPhylogeneticModel)

Computes specialized Fourier transform from the matrix that represents the full Fourier transform. 
This matrix transforms the reduced probability coordinates (corresponding to the equivalent classes) to the reduced Fourier coordinates.
"""
function fourier_transform(PM::GroupBasedPhylogeneticModel)
  FRp, p = Oscar.full_model_ring(phylogenetic_model(PM))
  FRq, q = Oscar.full_model_ring(PM)

  Rp, rp = Oscar.model_ring(phylogenetic_model(PM))
  Rq, rq = Oscar.model_ring(PM)

  p_classes = Oscar.equivalent_classes(phylogenetic_model(PM))
  q_classes = Oscar.equivalent_classes(PM)

  np = length(p_classes)
  nq = length(q_classes)
  
  ## We need to sort the equivalence classes: both inside each class as well as the collection of classes. 
  # keys_p_classes = collect(keys(p_classes))
  # sort!(keys_p_classes, rev = true)
  # keys_q_classes = collect(keys(q_classes))
  # sort!(keys_q_classes, rev = true)
  H = Rp.(hadamard(matrix_space(ZZ, n_states(PM), n_states(PM))))
  M = Rp.(Int.(zeros(nq, np)))
  
  for (i, q_index) in enumerate(sort(collect(keys(gens(Rq))); rev=true))
    current_fourier_class = q_classes[q_index]
    for (j, p_index) in enumerate(sort(collect(keys(gens(Rp))); rev=true))
      current_prob_class = p_classes[p_index]
      current_entries_in_M = [prod([H[y,x] for (x,y) in zip(FRp[pp], FRq[qq])]) for pp in current_prob_class, qq in current_fourier_class]
      M[i,j] = Rp.(1//(length(current_prob_class)*length(current_fourier_class))*sum(current_entries_in_M))
    end
  end
  
  
  return M
end

@doc raw"""
    coordinate_change(PM::GroupBasedPhylogeneticModel)

Constructs the map from probability coordinates to Fourier coordinates. This is a homomorphism between the model rings.
"""
function coordinate_change(PM::GroupBasedPhylogeneticModel)
  Rp, _ = Oscar.model_ring(phylogenetic_model(PM))
  Rq, _ = Oscar.model_ring(PM)
  
  M = fourier_transform(PM)
  hom(Rq, Rp, M * gens(_ring(Rp)))
end

@doc raw"""
    inverse_fourier_transform(PM::GroupBasedPhylogeneticModel)

Computes the matrix that transforms Fourier coordinates back to probability coordinates.
"""
function inverse_fourier_transform(PM::GroupBasedPhylogeneticModel)
  FRp, p = Oscar.full_model_ring(phylogenetic_model(PM))
  FRq, q = Oscar.full_model_ring(PM)

  Rp, _ = Oscar.model_ring(phylogenetic_model(PM))
  Rq, _ = Oscar.model_ring(PM)

  p_classes = Oscar.equivalent_classes(phylogenetic_model(PM))
  q_classes = Oscar.equivalent_classes(PM)

  np = length(p_classes)
  nq = length(q_classes)
  
  ## We need to sort the equivalence classes: both inside each class as well as the collection of classes. 
  # keys_p_classes = collect(keys(p_classes))
  # sort!(keys_p_classes, rev = true)
  # keys_q_classes = collect(keys(q_classes))
  # sort!(keys_q_classes, rev = true)
  H = Rq.(hadamard(matrix_space(ZZ, n_states(PM), n_states(PM))))
  Hinv = 1//n_states(PM) * H

  M = Rq.(Int.(zeros(np, nq)))
  for (i, p_index) in enumerate(sort(collect(keys(gens(Rp))); rev=true))
    current_prob_class = p_classes[p_index]
    for (j, q_index) in enumerate(sort(collect(keys(gens(Rq))); rev=true))
      current_fourier_class = q_classes[q_index]
      current_entries_in_M = [
        prod([Hinv[x,y] for (x,y) in zip(FRp[pp],FRq[qq])
                ]) for pp in current_prob_class, qq in current_fourier_class]
      M[i,j] = Rq.(sum(current_entries_in_M))
    end
  end
  return M
end

@doc raw"""
    inverse_coordinate_change(PM::GroupBasedPhylogeneticModel)

Constructs the map from Fourier coordinates to probability coordinates. This is the inverse of `coordinate_change`.
"""
function inverse_coordinate_change(PM::GroupBasedPhylogeneticModel)
  Rp, _ = Oscar.model_ring(phylogenetic_model(PM))
  Rq, _ = Oscar.model_ring(PM)
  
  M = inverse_fourier_transform(PM)
  hom(Rp, Rq, M*gens(_ring(Rq)))
end


###################################################################################
#
#       Construction of specific models
#
###################################################################################

@doc raw"""
    cavender_farris_neyman_model(G::AbstractGraph{Directed})
    cavender_farris_neyman_model(F::Field, G::AbstractGraph{Directed})

Constructs a `GroupBasedPhylogeneticModel` corresponding to the Cavender-Farris-Neyman model,
a 2-state model with a single transition probability per edge. 
If the field `F` is not provided, the model is constructed over the default rational field (`QQ`).

Example
```jldoctest
julia> cavender_farris_neyman_model(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]))
Group-based phylogenetic model on a tree with 3 leaves and 3 edges 
with root distribution [1//2, 1//2], 
transition matrices of the form 
 [:a :b;
  :b :a]
and fourier parameters of the form [:x, :y].
```

"""
function cavender_farris_neyman_model(F::Field, G::AbstractGraph{Directed})
  M = [:a :b;
       :b :a]
  x = [:x, :y] 

  group = collect(abelian_group(2))
  group = [group[1],group[2]]
  
  #PhylogeneticModel(G, M)
  GroupBasedPhylogeneticModel(F, G, M, x, group) 
end
cavender_farris_neyman_model(G::AbstractGraph{Directed}) = cavender_farris_neyman_model(QQ, G::AbstractGraph{Directed})

@doc raw"""
    jukes_cantor_model(G::AbstractGraph{Directed})
    jukes_cantor_model(F::Field, G::AbstractGraph{Directed})

Constructs a `GroupBasedPhylogeneticModel` corresponding to the Jukes-Cantor model,
a 4-state model where all non-diagonal transition probabilities are equal.
If the field `F` is not provided, the model is constructed over the default rational field (`QQ`).

Example
```jldoctest
julia> tree = phylogenetic_tree(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> jukes_cantor_model(tree)
Group-based phylogenetic model on a tree with 3 leaves and 3 edges 
with root distribution [1//4, 1//4, 1//4, 1//4], 
transition matrices of the form 
 [:a :b :b :b;
  :b :a :b :b;
  :b :b :a :b;
  :b :b :b :a]
and fourier parameters of the form [:x, :y, :y, :y].
```
"""
function jukes_cantor_model(F::Field, G::AbstractGraph{Directed})
  M = [:a :b :b :b;
       :b :a :b :b;
       :b :b :a :b;
       :b :b :b :a]
  x = [:x, :y, :y, :y] 
  GroupBasedPhylogeneticModel(F, G, M, x) 
end

jukes_cantor_model(G::AbstractGraph{Directed}) = jukes_cantor_model(QQ, G::AbstractGraph{Directed})

@doc raw"""
    kimura2_model(G::AbstractGraph{Directed})
    kimura2_model(F::Field, G::AbstractGraph{Directed})

Constructs a `GroupBasedPhylogeneticModel` corresponding to the Kimura 2-parameter model,
a 4-state model with two non-diagonal transition probabilities.
If the field `F` is not provided, the model is constructed over the default rational field (`QQ`).

Example
```jldoctest
julia> tree = phylogenetic_tree(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> kimura2_model(tree)
Group-based phylogenetic model on a tree with 3 leaves and 3 edges 
with root distribution [1//4, 1//4, 1//4, 1//4], 
transition matrices of the form 
 [:a :b :c :b;
  :b :a :b :c;
  :c :b :a :b;
  :b :c :b :a]
and fourier parameters of the form [:x, :y, :z, :z].
```
"""
function kimura2_model(F::Field, G::AbstractGraph{Directed})
  M = [:a :b :c :b;
       :b :a :b :c;
       :c :b :a :b;
       :b :c :b :a]
  x = [:x, :y, :z, :z] 
  
  GroupBasedPhylogeneticModel(F, G, M, x)
end

kimura2_model(G::AbstractGraph{Directed}) = kimura2_model(QQ, G::AbstractGraph{Directed})

@doc raw"""
    kimura3_model(G::AbstractGraph{Directed})
    kimura3_model(F::Field, G::AbstractGraph{Directed})

Constructs a `GroupBasedPhylogeneticModel` corresponding to the Kimura 3-parameter model,
a 4-state model with three non-diagonal different transition probabilities.
If the field `F` is not provided, the model is constructed over the default rational field (`QQ`).

Example
```jldoctest
julia> tree = phylogenetic_tree(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> kimura3_model(tree)
Group-based phylogenetic model on a tree with 3 leaves and 3 edges 
with root distribution [1//4, 1//4, 1//4, 1//4], 
transition matrices of the form 
 [:a :b :c :d;
  :b :a :d :c;
  :c :d :a :b;
  :d :c :b :a]
and fourier parameters of the form [:x, :y, :z, :t].
```
"""
function kimura3_model(F::Field, G::AbstractGraph{Directed})
  M = [:a :b :c :d;
       :b :a :d :c;
       :c :d :a :b;
       :d :c :b :a]
  x = [:x, :y, :z, :t] 
  
  GroupBasedPhylogeneticModel(F, G, M, x)
end

kimura3_model(G::AbstractGraph{Directed}) = kimura3_model(QQ, G::AbstractGraph{Directed})

@doc raw"""
    general_markov_model(G::AbstractGraph{Directed})
    general_markov_model(F::Field, G::AbstractGraph{Directed})

Constructs a `PhylogeneticModel` corresponding to a general Markov model with four states.
All transition matrix entries and root distribution entries are treated as independent parameters.
If the field `F` is not provided, the model is constructed over the default rational field (`QQ`).

Example
```jldoctest
julia> tree = phylogenetic_tree(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> general_markov_model(tree)
Phylogenetic model on a tree with 3 leaves and 3 edges
with root distribution [π1, π2, π3, π4] and transition matrices of the form
 [:m11 :m12 :m13 :m14;
  :m21 :m22 :m23 :m24;
  :m31 :m32 :m33 :m34;
  :m41 :m42 :m43 :m44].

```
"""
function general_markov_model(F::Field, G::AbstractGraph{Directed})

  M = [:m11 :m12 :m13 :m14;
       :m21 :m22 :m23 :m24;
       :m31 :m32 :m33 :m34;
       :m41 :m42 :m43 :m44]
  
  root_distr = [:π1, :π2, :π3, :π4]
  
  PhylogeneticModel(F, G, M, root_distr)
end

general_markov_model(G::AbstractGraph{Directed}) = general_markov_model(QQ, G::AbstractGraph{Directed})


@doc raw"""
    general_time_reversible_model(G::AbstractGraph{Directed})
    general_time_reversible_model(F::Field, G::AbstractGraph{Directed})

Constructs a `PhylogeneticModel` corresponding to a general time reversible model with four states.
All transition matrix entries and root distribution entries are treated as independent parameters.
If the field `F` is not provided, the model is constructed over the default rational field (`QQ`).

Example
```jldoctest
julia> tree = phylogenetic_tree(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> general_time_reversible_model(tree)
Phylogenetic model on a tree with 3 leaves and 3 edges
with root distribution [r[1], r[2], r[3], r[4]] and transition matrices of the form
 [r[1]*a r[2]*b r[3]*c r[4]*d;
  r[1]*b r[2]*a r[3]*e r[4]*f;
  r[1]*c r[2]*e r[3]*a r[4]*g;
  r[1]*d r[2]*f r[3]*g r[4]*a].
```
"""
function general_time_reversible_model(F::Field, G::AbstractGraph{Directed})
  Rp, r = polynomial_ring(QQ, :r => 1:4);
  RM, (a, b, c, d, e, f, g) = polynomial_ring(Rp, [:a, :b, :c, :d, :e, :f, :g]);

  M = [a*r[1] b*r[2] c*r[3] d*r[4];
       b*r[1] a*r[2] e*r[3] f*r[4];
       c*r[1] e*r[2] a*r[3] g*r[4];
       d*r[1] f*r[2] g*r[3] a*r[4]
  ];

  PM = PhylogeneticModel(F, G, M, r)
end

general_time_reversible_model(G::AbstractGraph{Directed}) = general_time_reversible_model(QQ, G::AbstractGraph{Directed})


