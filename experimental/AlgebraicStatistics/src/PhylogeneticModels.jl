################################################
#### PHYLOGENETIC DATA STRUCTURES & METHODS ####
################################################


###################################################################################
#
#       Phylogenetic Models
#
###################################################################################

@doc raw"""
    PhylogeneticModel{GT, M, L, T}
    PhylogeneticModel(F::Field, G::Graph{Directed}, trans_matrix_structure::Matrix, root_distribution::Union{Nothing, Vector} = nothing, varnames::VarName="p")

A data structure representing a general phylogenetic model on a directed tree.

This model defines a probability distribution on the states of the leaves of a phylogenetic tree,
based on a root distribution and transition matrices on the edges. The model is defined over a
base field (e.g., `QQ` for rational numbers).

**Type Parameters:**
* `GT`: The graph type, expected to be `Graph{Directed}` for a phylogenetic tree.
* `M`: The type of elements in the transition matrix structure, typically a `VarName` (a `Symbol`) or a polynomial ring element (`MPolyRingElem`).
* `L`: The type of the graph labelings, a `NamedTuple`.
* `T`: The type of elements in the root distribution vector.

**Fields:**
* `base_field::Field`: The base field of the model's parameters.
* `graph::GT`: The underlying directed graph (tree) of the model.
* `labelings::L`: A named tuple of graph maps for internal use.
* `trans_matrix_structure::Matrix{M}`: A matrix representing the symbolic structure of the transition matrices for each edge.
* `root_distribution::Vector{T}`: A vector representing the distribution at the root of the tree.
* `n_states::Int`: The number of states in the model.
* `model_parameter_name::VarName`: The symbolic name for the variables of the model ring.
"""
@attributes mutable struct PhylogeneticModel{GT, M, L, T} <: GraphicalModel{GT, L}
  # do need to for T to be directed here? YES!
  base_field::Field
  graph::GT
  labelings::L
  trans_matrix_structure::Matrix{M}
  root_distribution::Vector{T}
  n_states::Int
  model_parameter_name::VarName

  function PhylogeneticModel(F::Field,
                             G::Graph{Directed},
                             trans_matrix_structure::Matrix,
                             root_distribution::Union{Nothing, Vector} = nothing,
                             varnames::VarName="p")
    n_states = size(trans_matrix_structure)[1]
    if isnothing(root_distribution)
      root_distribution = F.(repeat([1//n_states], outer = n_states))
    end
    graph_maps = NamedTuple(_graph_maps(G))
    graph_maps = isempty(graph_maps) ? nothing : graph_maps
    return new{Graph{Directed}, typeof(first(trans_matrix_structure)), 
               typeof(graph_maps), typeof(first(root_distribution))}(F, G,
                                                                    graph_maps,
                                                                    trans_matrix_structure,
                                                                    root_distribution,
                                                                    n_states,
                                                                    varnames)
  end

  function PhylogeneticModel(G::Graph{Directed},
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

  function PhylogeneticModel(G::Graph{Directed},
                             trans_matrix_structure::Matrix,
                             root_distribution::Union{Nothing, Vector} = nothing,
                             varnames::VarName="p")
    return PhylogeneticModel(QQ, G, trans_matrix_structure, 
                             root_distribution, varnames)
  end
end

@doc raw"""
      phylogenetic_model(F::Field, G::Graph{Directed}, trans_matrix_structure::Matrix, root_distribution::Union{Nothing, Vector} = nothing, varnames::VarName="p")
      phylogenetic_model(G::Graph{Directed}, trans_matrix_structure::Matrix, root_distribution::Union{Nothing, Vector} = nothing, varnames::VarName="p")
      phylogenetic_model(G::Graph{Directed}, trans_matrix_structure::Matrix{M}, root_distribution::Vector{T}, varnames::VarName="p") where {M <: MPolyRing{<:MPolyRingElem}, T <: MPolyRing}
      
Construct a `PhylogeneticModel` from a field `F`, a directed graph `G`, a symbolic transition matrix `trans_matrix_structure`, and an optional `root_distribution`.
If `root_distribution` is not provided, a uniform distribution is assumed. If `F` is not provided, it constructs a `PhylogeneticModel` using the default rational field (`QQ`).

If `trans_matrix_structure` is a `Matrix{<: MPolyRing{<:MPolyRingElem}}` and `root_distribution` is a `Vector{<: MPolyRing}`, then the constructor ensures the rings of the root distribution parameters and the transition matrices are compatible.
"""
function phylogenetic_model(F::Field, G::Graph{Directed},
                            trans_matrix_structure::Matrix,
                            root_distribution::Union{Nothing, Vector} = nothing,
                            varnames::VarName="p")
  PhylogeneticModel(F, G, trans_matrix_structure, root_distribution, varnames)
end

function phylogenetic_model(G::Graph{Directed}, trans_matrix_structure::Matrix{M},
                            root_distribution::Vector{T},
                            varnames::VarName="p") where {M <: MPolyRing{<:MPolyRingElem}, T <: MPolyRing}
  PhylogeneticModel(G, trans_matrix_structure, root_distribution, varnames)
end

function phylogenetic_model(G::Graph{Directed}, trans_matrix_structure::Matrix,
                            root_distribution::Union{Nothing, Vector} = nothing,
                            varnames::VarName="p")
  return PhylogeneticModel(G, trans_matrix_structure, root_distribution, varnames)
end


@doc raw"""
    n_states(PM::PhylogeneticModel)

Return the number of states in the phylogenetic model.

# Example
```jldoctest
julia> PM = general_markov_model(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> n_states(PM)
4
```
"""
n_states(PM::PhylogeneticModel) = PM.n_states

@doc raw"""
    transition_matrix(PM::PhylogeneticModel)

Return the structure of the transition matrices of the model.

# Example
```jldoctest
julia> PM = general_markov_model(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> transition_matrix(PM)
4×4 Matrix{Symbol}:
 :m11  :m12  :m13  :m14
 :m21  :m22  :m23  :m24
 :m31  :m32  :m33  :m34
 :m41  :m42  :m43  :m44

```
"""
transition_matrix(PM::PhylogeneticModel) = PM.trans_matrix_structure

@doc raw"""
    root_distribution(PM::PhylogeneticModel)

Return the root distribution vector of the model.

# Example
```jldoctest
julia> PM = general_markov_model(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> root_distribution(PM)
4-element Vector{Symbol}:
 :π1
 :π2
 :π3
 :π4
```
"""
root_distribution(PM::PhylogeneticModel) = PM.root_distribution

@doc raw"""
    base_field(PM::PhylogeneticModel)

Return the base field of the model's rings.

# Example
```jldoctest
julia> PM = general_markov_model(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> base_field(PM)
Rational field
```
"""
base_field(PM::PhylogeneticModel) = PM.base_field

@doc raw"""
    varnames(PM::PhylogeneticModel)

Return the symbolic name for the probability coordinates.

# Example
```jldoctest
julia> PM = general_markov_model(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> varnames(PM)
"p"
```
"""
varnames(PM::PhylogeneticModel) = PM.model_parameter_name


@doc raw"""
    parameter_ring(PM::PhylogeneticModel{Graph{Directed}, <: VarName, L, T}; cached=false) where {L, T <: FieldElem}
    parameter_ring(PM::PhylogeneticModel{Graph{Directed}, <: VarName}; cached=false)
    parameter_ring(PM::PhylogeneticModel{Graph{Directed}, <: MPolyRingElem, L, T}; cached=false) where {L, T <: FieldElem}
    parameter_ring(PM::PhylogeneticModel{Graph{Directed}, <: MPolyRingElem, L, T}; cached=false) where {L, T <: MPolyRingElem}
    parameter_ring(PM::PhylogeneticModel{Graph{Directed}, <: MPolyRingElem, L, T}; cached=false) where {L, T <: AbstractAlgebra.Generic.RationalFunctionFieldElem}

Create the polynomial ring for the parameters of a phylogenetic model `PhylogeneticModel{GT, M, L, T}`.

If `M  <: VarName`, returns a tuple containing:
1.  The polynomial ring.
2.  A dictionary mapping parameter variables and edges to the corresponding ring generators.
3.  The root distribution vector.

If `M  <: MPolyRingElem`, then returns a tuple containing:
1.  The polynomial ring.
2.  A dictionary mapping each edge to a homomorphism (`Oscar.MPolyAnyMap`) from the original
    transition matrix ring to the new parameter ring.
3.  The root distribution vector.
"""
@attr Tuple{
  MPolyRing, 
  GraphTransDict,
  Vector{T}
} function parameter_ring(PM::PhylogeneticModel{Graph{Directed}, <: VarName, L, T}; cached=false) where {L, T <: FieldElem} 
  vars = unique(transition_matrix(PM))
  edge_gens = [x => 1:n_edges(graph(PM)) for x in vars]
  R, x... = polynomial_ring(base_field(PM), edge_gens...; cached=cached)
  
  R, Dict{Tuple{VarName, Edge}, MPolyRingElem}(
    (vars[i], e) => x[i][j] for i in 1:length(vars), 
      (j,e) in enumerate(sort_edges(graph(PM)))
      ), root_distribution(PM)
end

@attr Tuple{
  MPolyRing,
  GraphTransDict,
  Vector{<:MPolyRingElem}
} function parameter_ring(PM::PhylogeneticModel{Graph{Directed}, <: VarName}; cached=false)
  vars = unique(transition_matrix(PM))
  edge_gens = [x => 1:n_edges(graph(PM)) for x in vars]
  R, r, x... = polynomial_ring(base_field(PM),
                               root_distribution(PM),
                               edge_gens...; cached=cached)

  R, Dict{Tuple{VarName, Edge}, MPolyRingElem}(
    (vars[i], e) => x[i][j] for i in 1:length(vars), 
      (j,e) in enumerate(sort_edges(graph(PM)))
      ), r
end

@attr Tuple{
  MPolyRing, 
  Dict{Edge, Oscar.MPolyAnyMap}, 
  Vector{T}
} function parameter_ring(PM::PhylogeneticModel{Graph{Directed}, <: MPolyRingElem, L, T}; cached=false) where {L, T <: FieldElem}
  trans_ring = parent(first(transition_matrix(PM)))
  transition_vars = gens(trans_ring)

  edge_gens = [x => 1:n_edges(graph(PM)) for x in Symbol.(transition_vars)]
  R, x... = polynomial_ring(base_field(PM), edge_gens...; cached=cached)
  
  dict_maps = Dict{Edge, Oscar.MPolyAnyMap}()
  for (j,e) in enumerate(Oscar.sort_edges(graph(PM)))
      map = [x[i][j] for i in 1:length(transition_vars)]
      dict_maps[e] = hom(trans_ring, R, map)
  end

  R, dict_maps, root_distribution(PM)
end

@attr Tuple{
  MPolyRing, 
  Dict{Edge, Oscar.MPolyAnyMap}, 
  Vector{T}
} function parameter_ring(PM::PhylogeneticModel{Graph{Directed}, <: MPolyRingElem, L, T}; cached=false) where {L, T <: MPolyRingElem}
  trans_ring = parent(first(transition_matrix(PM)))
  transition_vars = gens(trans_ring)
  root_vars = gens(coefficient_ring(trans_ring))

  edge_gens = [x => 1:n_edges(graph(PM)) for x in Symbol.(transition_vars)]
  R, rv, x... = polynomial_ring(base_field(PM), Symbol.(root_vars), edge_gens..., ; cached=cached)

  coef_map = hom(coefficient_ring(trans_ring), R, rv)

  dict_maps = Dict{Edge, Oscar.MPolyAnyMap}()
  for (j,e) in enumerate(Oscar.sort_edges(graph(PM)))
      map = [x[i][j] for i in 1:length(transition_vars)]
      dict_maps[e] = hom(trans_ring, R, coef_map, map)
  end

  R, dict_maps, rv

end

@attr Tuple{
  MPolyRing{T}, 
  Dict{Edge, Oscar.MPolyAnyMap}, 
  Vector{T}
} function parameter_ring(PM::PhylogeneticModel{Graph{Directed}, <: MPolyRingElem, L, T}; cached=false) where {L, T <: AbstractAlgebra.Generic.RationalFunctionFieldElem}
  trans_ring = parent(first(transition_matrix(PM)))
  transition_vars = gens(trans_ring)
  root_vars = gens(coefficient_ring(trans_ring))

  edge_gens = [x => 1:n_edges(graph(PM)) for x in Symbol.(transition_vars)]
  R, x... = polynomial_ring(coefficient_ring(trans_ring), edge_gens..., ; cached=cached)

  dict_maps = Dict{Edge, Oscar.MPolyAnyMap}()
  for (j,e) in enumerate(Oscar.sort_edges(graph(PM)))
    map = [x[i][j] for i in 1:length(transition_vars)]
    dict_maps[e] = hom(trans_ring, R, map)
  end
  R, dict_maps, root_vars
end

function Base.show(io::IO, PM::PhylogeneticModel)
  gr = graph(PM)

  nl = length(leaves(gr))
  ne = length(collect(edges(gr)))
  root_dist = join(root_distribution(PM), ", " )
 
  print(io, "Phylogenetic model on a tree with $(nl) leaves and $(ne) edges \n") # \n )
  print(io, "with root distribution [$(root_dist)] \n")
  print(io, "and transition matrices of the form \n ")

  M = string(PM.trans_matrix_structure)
  M = split(M, "[", limit=2)[2]
  print(io, "[", replace(M, ";" => ";\n "))
  print(io, ". ")
end


###################################################################################
#
#       Group-based phylogenetic models
#
###################################################################################

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
  phylo_model::PhylogeneticModel{GT, <: VarName, L, <: RingElem} 
  fourier_param_structure::Vector{<: VarName}
  group::Vector{FinGenAbGroupElem}
  model_parameter_name::VarName
  
  function GroupBasedPhylogeneticModel(F::Field, 
                                       G::Graph{Directed},
                                       trans_matrix_structure::Matrix{<: VarName},
                                       fourier_param_structure::Vector{<: VarName},
                                       group::Union{Nothing, Vector{FinGenAbGroupElem}} = nothing,
                                       root_distribution::Union{Nothing, Vector} = nothing,
                                       varnames_phylo_model::VarName="p",
                                       varnames_group_based::VarName="q")
    if isnothing(group)
      group = collect(abelian_group(2,2))
      group = [group[1],group[3],group[2],group[4]]
    end

    graph_maps = NamedTuple(_graph_maps(G))
    graph_maps = isempty(graph_maps) ? nothing : graph_maps
    return new{typeof(G), typeof(graph_maps)}(PhylogeneticModel(F, G, trans_matrix_structure, 
                                                                root_distribution,
                                                                varnames_phylo_model),
                                              fourier_param_structure,
                                              group, varnames_group_based)
  end

  function GroupBasedPhylogeneticModel(G::Graph{Directed},
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
    group_based_phylogenetic_model(F::Field, G::Graph{Directed}, trans_matrix_structure::Matrix{<: VarName}, fourier_param_structure::Vector{<: VarName}, group::Union{Nothing, Vector{FinGenAbGroupElem}} = nothing, root_distribution::Union{Nothing, Vector} = nothing, varnames_phylo_model::VarName="p", varnames_group_based::VarName="q")
    group_based_phylogenetic_model(G::Graph{Directed}, trans_matrix_structure::Matrix{<: VarName}, fourier_param_structure::Vector{<: VarName}, group::Union{Nothing, Vector{FinGenAbGroupElem}} = nothing, root_distribution::Union{Nothing, Vector} = nothing, varnames_phylo_model::VarName="p", varnames_group_based::VarName="q")


Construct a `GroupBasedPhylogeneticModel` from a field `F`, a directed graph `G`, a symbolic transition matrix, symbolic Fourier parameters, an optional group, and optional root distribution. If `F` is not provided it constructs a `GroupBasedPhylogeneticModel` using the default rational field (`QQ`).
"""

function group_based_phylogenetic_model(F::Field, 
                                      G::Graph{Directed},
                                      trans_matrix_structure::Matrix{<: VarName},
                                      fourier_param_structure::Vector{<: VarName},
                                      group::Union{Nothing, Vector{FinGenAbGroupElem}} = nothing,
                                      root_distribution::Union{Nothing, Vector} = nothing,
                                      varnames_phylo_model::VarName="p",
                                      varnames_group_based::VarName="q")
   GroupBasedPhylogeneticModel(F, G, trans_matrix_structure, fourier_param_structure, group,
                               root_distribution, varnames_phylo_model, varnames_group_based)
end

function group_based_phylogenetic_model(G::Graph{Directed},
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

@doc raw"""
    graph(PM::GroupBasedPhylogeneticModel)

Return the graph of the underlying phylogenetic model.

# Example
```jldoctest
julia> PM = jukes_cantor_model(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> graph(PM)
Directed graph with 4 nodes and the following edges:
(4, 1)(4, 2)(4, 3)

```
"""
graph(PM::GroupBasedPhylogeneticModel) = PM.phylo_model.graph # ? Antony

@doc raw"""
    phylogenetic_model(PM::GroupBasedPhylogeneticModel)

Return the underlying `PhylogeneticModel`.

# Example
```jldoctest
julia> PM = jukes_cantor_model(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> phylogenetic_model(PM)
Phylogenetic model on a tree with 3 leaves and 3 edges 
with root distribution [1//4, 1//4, 1//4, 1//4] 
and transition matrices of the form 
 [:a :b :b :b;
  :b :a :b :b;
  :b :b :a :b;
  :b :b :b :a]. 
```
"""
phylogenetic_model(PM::GroupBasedPhylogeneticModel) = PM.phylo_model


@doc raw"""
    n_states(PM::GroupBasedPhylogeneticModel)

Return the number of states of the underlying phylogenetic model.

# Example
```jldoctest
julia> PM = jukes_cantor_model(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> n_states(PM)
4
```
"""
n_states(PM::GroupBasedPhylogeneticModel) = PM.phylo_model.n_states

@doc raw"""
    transition_matrix(PM::GroupBasedPhylogeneticModel)

Return the structure of the transition matrices of the model.

# Example
```jldoctest
julia> PM = jukes_cantor_model(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> transition_matrix(PM)
4×4 Matrix{Symbol}:
 :a  :b  :b  :b
 :b  :a  :b  :b
 :b  :b  :a  :b
 :b  :b  :b  :a
```
"""
transition_matrix(PM::GroupBasedPhylogeneticModel) = PM.phylo_model.trans_matrix_structure

@doc raw"""
    base_field(PM::GroupBasedPhylogeneticModel)

Return the base field of the underlying phylogenetic model.

# Example
```jldoctest
julia> PM = jukes_cantor_model(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> base_field(PM)
Rational field
```
"""
base_field(PM::GroupBasedPhylogeneticModel) = PM.phylo_model.base_field

@doc raw"""
    root_distribution(PM::GroupBasedPhylogeneticModel)

Return the root distribution of the underlying phylogenetic model.

# Example
```jldoctest
julia> PM = jukes_cantor_model(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> root_distribution(PM)
4-element Vector{QQFieldElem}:
 1//4
 1//4
 1//4
 1//4
```
"""
root_distribution(PM::GroupBasedPhylogeneticModel) = PM.phylo_model.root_distribution

@doc raw"""
    group(PM::GroupBasedPhylogeneticModel)

Return the finite abelian group associated with the model.

# Example
```jldoctest
julia> PM = jukes_cantor_model(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

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
julia> PM = jukes_cantor_model(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> fourier_parameters(PM)
4-element Vector{Symbol}:
 :x
 :y
 :y
 :y
```
"""
fourier_parameters(PM::GroupBasedPhylogeneticModel) = PM.fourier_param_structure

varnames_probabilities(PM::GroupBasedPhylogeneticModel) = PM.phylo_model.model_parameter_name # ? Antony

varnames_fourier(PM::GroupBasedPhylogeneticModel) = PM.model_parameter_name # ? Antony

@doc raw"""
    varnames(PM::GroupBasedPhylogeneticModel)

Return the variable name for the Fourier coordinates (alias for `varnames_fourier`).

# Example
```jldoctest
julia> PM = jukes_cantor_model(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]));

julia> varnames(PM)
"q"
```
"""
varnames(PM::GroupBasedPhylogeneticModel) = PM.model_parameter_name


@doc raw"""
    parameter_ring(PM::GroupBasedPhylogeneticModel; cached=false)

Create the polynomial ring for the Fourier parameters of the model.

Returns a tuple containing:
1.  The polynomial ring.
2.  A dictionary mapping parameter variables and edges to the corresponding ring generators.
"""
@attr Tuple{
  MPolyRing,
  GraphTransDict
} function parameter_ring(PM::GroupBasedPhylogeneticModel; cached=false)
  vars = unique(fourier_parameters(PM))
  edge_gens = [x => 1:n_edges(graph(PM)) for x in vars]
  R, x... = polynomial_ring(base_field(PM), edge_gens...; cached=cached)

  R, Dict{Tuple{VarName, Edge}, MPolyRingElem}(
    (vars[i], e) => x[i][j] for i in 1:length(vars), (j,e) in enumerate(sort_edges(graph(PM)))
  )
end

function Base.show(io::IO, PM::GroupBasedPhylogeneticModel)

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


###################################################################################
#
#       Parametrizations
#
###################################################################################

@doc raw"""
    full_model_ring(PM::Union{PhylogeneticModel, GroupBasedPhylogeneticModel}; cached=false)

Creates the full model ring for the probability distributions at the leaves (probability coordinates) for a `PhylogeneticModel`.
Equivalently, creates the full model ring in the Fourier coordinates for a `GroupBasedPhylogeneticModel`.
This ring has a generator for every possible state configuration at the leaves.
"""
@attr Tuple{
  ModelRing{T, U}, 
  Dict{T, U}
} where {T, U} function full_model_ring(PM::Union{PhylogeneticModel, GroupBasedPhylogeneticModel}; cached=false)
  l_indices = leaves_indices(PM)

  return model_ring(base_field(PM), varnames(PM) => l_indices; cached=cached)
end

@doc raw"""
    full_parametrization(PM::PhylogeneticModel)

Constructs the parametrization map from the full model ring in _probability_ coordinates (with generators for all leaf probability
configurations) to the parameter ring (with generators from the transition matrices and the root distribution).
"""
@attr MPolyAnyMap function full_parametrization(PM::PhylogeneticModel)
  gr = graph(PM)

  R, _ = full_model_ring(PM)
  S, _ = parameter_ring(PM)

  lvs_indices = leaves_indices(PM::PhylogeneticModel)

  map = [leaves_probability(PM, Dict(i => k[i] for i in 1:n_leaves(gr))) for k in lvs_indices]
  hom(R, S, reduce(vcat, map))

end

@doc raw"""
    full_parametrization(PM::GroupBasedPhylogeneticModel)

Constructs the parametrization map from the full model ring in _Fourier_ coordinates (with generators for all leaf probability
configurations) to the parameter ring (with generators from the Fourier parameters).
"""
function full_parametrization(PM::GroupBasedPhylogeneticModel)
  gr = graph(PM)

  R, _ = full_model_ring(PM)
  S, _ = parameter_ring(PM)

  lvs_indices = leaves_indices(PM)

  map = [leaves_fourier(PM, Dict(i => k[i] for i in 1:n_leaves(gr))) for k in lvs_indices]
  hom(R, S, reduce(vcat, map))

end

@doc raw"""
    equivalent_classes(PM::Union{PhylogeneticModel, GroupBasedPhylogeneticModel})

Calculates the equivalence classes of the model's parametrizations. Equivalent polynomials
are those that are identical under the chosen parametrization. This is a crucial step for
dimension reduction in the model's algebraic variety.
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
    parametrization(PM::PhylogeneticModel)

Constructs the parametrization map from the reduced model ring of probability coordinates (with generators for
each equivalence class) to the parameter ring.
"""
function parametrization(PM::PhylogeneticModel)
  R, _ = model_ring(PM)
  S, _ = parameter_ring(PM)

  # need to double check the order of the keys here is consitent,
  # otherwise we need to add a sort
  lvs_indices = keys(gens(R))
  map = [leaves_probability(PM, Dict(i => k[i] for i in 1:n_leaves(graph(PM)))) for k in lvs_indices]
  hom(R, S, reduce(vcat, map))
end

@doc raw"""
    parametrization(PM::GroupBasedPhylogeneticModel)

Constructs the parametrization map from the reduced model of Fourier coordinates ring to the Fourier parameter ring.
"""
function parametrization(PM::GroupBasedPhylogeneticModel)
  R, _ = model_ring(PM)
  S, _ = parameter_ring(PM)
  
  lvs_indices = keys(gens(R))

  map = [leaves_fourier(PM, Dict(i => k[i] for i in 1:n_leaves(graph(PM)))) for k in lvs_indices]
  hom(R, S, reduce(vcat, map))
end

@doc raw"""
    affine_parametrization(PM::PhylogeneticModel)

Constructs an affine parametrization of the model by imposing the transition matrices rows sum to one and the root distribution sums to one.
"""
function affine_parametrization(PM::PhylogeneticModel)
  R, p = model_ring(PM) 
  S, _ = parameter_ring(PM)
  f = parametrization(PM)
  ind = keys(gens(R))

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
"""
function affine_parametrization(PM::GroupBasedPhylogeneticModel)
  R, q = model_ring(PM)
  S, _ = parameter_ring(PM)
  f = parametrization(PM)
  ind = keys(gens(R))

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
    cavender_farris_neyman_model(G::Graph{Directed})

Constructs a `GroupBasedPhylogeneticModel` corresponding to the Cavender-Farris-Neyman model,
a 2-state model with a single transition probability per edge.

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
function cavender_farris_neyman_model(G::Graph{Directed})
  M = [:a :b;
       :b :a]
  x = [:x, :y] 

  group = collect(abelian_group(2))
  group = [group[1],group[2]]
  
  #PhylogeneticModel(G, M)
  GroupBasedPhylogeneticModel(G, M, x, group) 
end

@doc raw"""
    jukes_cantor_model(G::Graph{Directed})

Constructs a `GroupBasedPhylogeneticModel` corresponding to the Jukes-Cantor model,
a 4-state model where all non-diagonal transition probabilities are equal.

Example
```jldoctest
julia> jukes_cantor_model(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]))
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
function jukes_cantor_model(G::Graph{Directed})
  M = [:a :b :b :b;
       :b :a :b :b;
       :b :b :a :b;
       :b :b :b :a]
  # x = [Symbol("x[1]"), Symbol("x[2]"), Symbol("x[2]"), Symbol("x[2]")]
  # x = [:x1, :x2, :x2, :x2] #? Antony
  x = [:x, :y, :y, :y] 
  
  #PhylogeneticModel(G, M)
  GroupBasedPhylogeneticModel(G, M, x) 
end

@doc raw"""
    kimura2_model(G::Graph{Directed})

Constructs a `GroupBasedPhylogeneticModel` corresponding to the Kimura 2-parameter model,
a 4-state model with two non-diagonal transition probabilities.

Example
```jldoctest
julia> kimura2_model(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]))
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
function kimura2_model(G::Graph{Directed})
  M = [:a :b :c :b;
       :b :a :b :c;
       :c :b :a :b;
       :b :c :b :a]
  x = [:x, :y, :z, :z] 
  
  GroupBasedPhylogeneticModel(G, M, x)
end

@doc raw"""
    kimura3_model(G::Graph{Directed})

Constructs a `GroupBasedPhylogeneticModel` corresponding to the Kimura 3-parameter model,
a 4-state model with three non-diagonal different transition probabilities.

Example
```jldoctest
julia> kimura3_model(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]))
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
function kimura3_model(G::Graph{Directed})
  M = [:a :b :c :d;
       :b :a :d :c;
       :c :d :a :b;
       :d :c :b :a]
  x = [:x, :y, :z, :t] 
  
  GroupBasedPhylogeneticModel(G, M, x)
end
  
@doc raw"""
    general_markov_model(G::Graph{Directed})

Constructs a `PhylogeneticModel` corresponding to a general Markov model with four states.
All transition matrix entries and root distribution entries are treated as independent parameters.

Example
```jldoctest
julia> general_markov_model(graph_from_edges(Directed,[[4,1], [4,2], [4,3]]))
Phylogenetic model on a tree with 3 leaves and 3 edges 
with root distribution [π1, π2, π3, π4] 
and transition matrices of the form 
 [:m11 :m12 :m13 :m14;
  :m21 :m22 :m23 :m24;
  :m31 :m32 :m33 :m34;
  :m41 :m42 :m43 :m44]. 
```
"""
function general_markov_model(G::Graph{Directed})

  M = [:m11 :m12 :m13 :m14;
       :m21 :m22 :m23 :m24;
       :m31 :m32 :m33 :m34;
       :m41 :m42 :m43 :m44]
  
  root_distr = [:π1, :π2, :π3, :π4]
  
  PhylogeneticModel(G, M, root_distr)
end

@doc raw"""
    general_markov_model(G::Graph{Directed})

Constructs a `PhylogeneticModel` corresponding to a general Markov model with four states.
All transition matrix entries and root distribution entries are treated as independent parameters.

Example
```jldoctest
julia> general_time_reversible_model(tree)
Phylogenetic model on a tree with 3 leaves and 3 edges 
with root distribution [r[1], r[2], r[3], r[4]] 
and transition matrices of the form 
 [r[1]*a r[2]*b r[3]*c r[4]*d;
  r[1]*b r[2]*a r[3]*e r[4]*f;
  r[1]*c r[2]*e r[3]*a r[4]*g;
  r[1]*d r[2]*f r[3]*g r[4]*a]. 
```
"""
function general_time_reversible_model(G::Graph{Directed})

  Rp, r = polynomial_ring(QQ, :r => 1:4);
  RM, (a, b, c, d, e, f, g) = polynomial_ring(Rp, [:a, :b, :c, :d, :e, :f, :g]);

  M = [a*r[1] b*r[2] c*r[3] d*r[4];
       b*r[1] a*r[2] e*r[3] f*r[4];
       c*r[1] e*r[2] a*r[3] g*r[4];
       d*r[1] f*r[2] g*r[3] a*r[4]
  ];

  PM = PhylogeneticModel(G, M, r)
end

