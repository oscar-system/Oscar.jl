```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Algebraic Phylogenetics

The Algebraic Phylogenetics module in OSCAR provides a comprehensive toolkit for the algebraic and geometric study of phylogenetic models. It allows for the construction of various evolutionary models on phylogenetic trees, the creation of their associated polynomial rings and parametrizations, and the analysis of their algebraic properties, such as phylogenetic invariants.

## Model Constructors

Two main types of models can be constructed: a general `PhylogeneticModel` and a more specialized `GroupBasedPhylogeneticModel`.

#### `PhylogeneticModel`
A `PhylogeneticModel` is defined by a directed tree, a base field, the symbolic structure of the transition matrices for each edge, and a probability distribution at the root.
```@docs
PhylogeneticModel(F::Field, G::Graph{Directed}, trans_matrix_structure::Matrix, root_distribution::Union{Nothing, Vector} = nothing, varname::VarName="p")
PhylogeneticModel(G::Graph{Directed},trans_matrix_structure::Matrix{M},root_distribution::Vector{T},varname::VarName="p") where {M <:  MPolyRing{<:MPolyRingElem}, T <: MPolyRing}
PhylogeneticModel(G::Graph{Directed}, trans_matrix_structure::Matrix, root_distribution::Union{Nothing, Vector} = nothing, varname::VarName="p")
```

Here are some examples of how to construct phylogenetic models:
```jldoctest Example_F81:
julia> tree = graph_from_edges(Directed,[[4,1], [4,2], [4,3]]);


# Jukes-Cantor model with symbolic parameters
julia> M_JC = [:a :b :b :b;
               :b :a :b :b;
               :b :b :a :b;
               :b :b :b :a];

julia> PM_jukesCantor = PhylogeneticModel(tree, M_JC)
Phylogenetic model on a tree with 3 leaves and 3 edges 
with root distribution [1//4, 1//4, 1//4, 1//4] 
and transition matrices of the form 
 [:a :b :b :b;
  :b :a :b :b;
  :b :b :a :b;
  :b :b :b :a]. 


# Felsenstein 1981 model with polynomial parameters
julia> Rp, r = rational_function_field(QQ, :r => 1:4);

julia> RM, (a, b) = polynomial_ring(Rp, [:a, :b]);

julia> M_F81 = [a*r[1] b*r[2] b*r[3] b*r[4];
                b*r[1] a*r[2] b*r[3] b*r[4];
                b*r[1] b*r[2] a*r[3] b*r[4];
                b*r[1] b*r[2] b*r[3] a*r[4]];

julia> PM_Felsenstein81 = PhylogeneticModel(tree, M_F81, r)
Phylogenetic model on a tree with 3 leaves and 3 edges 
with root distribution [r[1], r[2], r[3], r[4]] 
and transition matrices of the form 
 [r[1]*a r[2]*b r[3]*b r[4]*b;
  r[1]*b r[2]*a r[3]*b r[4]*b;
  r[1]*b r[2]*b r[3]*a r[4]*b;
  r[1]*b r[2]*b r[3]*b r[4]*a]. 
```

#### `GroupBasedPhylogeneticModel`

For models exhibiting symmetries that can be described by a finite abelian group, the `GroupBasedPhylogeneticModel` is used. This structure requires the symbolic Fourier parameters and the associated group, in addition to the standard phylogenetic model components.

```@docs
GroupBasedPhylogeneticModel(F::Field, G::Graph{Directed},trans_matrix_structure::Matrix{<: VarName},fourier_param_structure::Vector{<: VarName},group::Union{Nothing, Vector{FinGenAbGroupElem}} = nothing,root_distribution::Union{Nothing, Vector} = nothing,varname_phylo_model::VarName="p",varname_group_based::VarName="q")
GroupBasedPhylogeneticModel(G::Graph{Directed},trans_matrix_structure::Matrix{<: VarName},fourier_param_structure::Vector{<: VarName},group::Union{Nothing, Vector{FinGenAbGroupElem}} = nothing,root_distribution::Union{Nothing, Vector} = nothing,varname_phylo_model::VarName="p",varname_group_based::VarName="q")
```

For example, the Jukes-Cantor model can be defined as a group-based model:
```jldoctest Example_JC:
julia> tree = graph_from_edges(Directed,[[4,1], [4,2], [4,3]]);

julia> M = [:a :b :b :b;
            :b :a :b :b;
            :b :b :a :b;
            :b :b :b :a];

julia>   x = [:x, :y, :y, :y]; 

julia>   GroupBasedPhylogeneticModel(tree, M, x) 
Phylogenetic model on a tree with 3 leaves and 3 edges 
with root distribution [1//4, 1//4, 1//4, 1//4], 
transition matrices of the form 
 [:a :b :b :b;
  :b :a :b :b;
  :b :b :a :b;
  :b :b :b :a]
and fourier parameters of the form [:x, :y, :y, :y].
```


### Classic phylogenetic models
Several classic phylogenetic models are pre-defined for convenience.
```@docs
cavender_farris_neyman_model(graph::Graph{Directed})
jukes_cantor_model(graph::Graph{Directed})
kimura2_model(graph::Graph{Directed})
kimura3_model(graph::Graph{Directed})
general_markov_model(graph::Graph{Directed})
```

## Properties
You can access the properties of a model using the following accessor functions.

#### `PhylogeneticModel` Properties

```@docs
n_states(PM::PhylogeneticModel)
transition_matrix(PM::PhylogeneticModel)
root_distribution(PM::PhylogeneticModel)
base_field(PM::PhylogeneticModel)
varname(PM::PhylogeneticModel)
```

#### `GroupBasedPhylogeneticModel` Properties

Objects of type `GroupBasedPhylogeneticModel` are built upon a `PhylogeneticModel`, which can be accessed directly:
```@docs
phylogenetic_model(PM::GroupBasedPhylogeneticModel)
```

Additionally, `GroupBasedPhylogeneticModel` has its own properties, many of which are aliases for the underlying model's properties.

```@docs
n_states(PM::GroupBasedPhylogeneticModel)
transition_matrix(PM::GroupBasedPhylogeneticModel)
root_distribution(PM::GroupBasedPhylogeneticModel) 
base_field(PM::GroupBasedPhylogeneticModel)
group(PM::GroupBasedPhylogeneticModel)
fourier_parameters(PM::GroupBasedPhylogeneticModel)
varname(PM::GroupBasedPhylogeneticModel)
```

## Phylogenetic rings
Two main rings are associated with any phylogenetic model: the _parameter ring_, which contains the unobserved parameters of the evolutionary process (e.g., transition probabilities), and the _model ring_, which represents the observed data (e.g., joint probability distribution at the leaves).

### Parameter ring
The parameter ring's generators correspond to the symbolic parameters of the model. For a `PhylogeneticModel`, these are the entries of the transition matrices and root distribution. For a `GroupBasedPhylogeneticModel`, they are the Fourier parameters.

```@docs
parameter_ring(PM::PhylogeneticModel{Graph{Directed}, <: VarName, L, T}; cached=false) where {L, T}
parameter_ring(PM::PhylogeneticModel{Graph{Directed}, <: MPolyRingElem, L, T}; cached=false) where {L, T}
parameter_ring(PM::GroupBasedPhylogeneticModel; cached=false)
```

For example

```jldoctest Example_JC:
julia> tree = graph_from_edges(Directed,[[4,1], [4,2], [4,3]]);

julia> pm = kimura2_model(tree)
Phylogenetic model on a tree with 3 leaves and 3 edges 
with root distribution [1//4, 1//4, 1//4, 1//4], 
transition matrices of the form 
 [:a :b :c :b;
  :b :a :b :c;
  :c :b :a :b;
  :b :c :b :a]
and fourier parameters of the form [:x, :y, :z, :z].

julia> S, S_gens = parameter_ring(pm);

julia> S
Multivariate polynomial ring in 9 variables x[1], x[2], x[3], y[1], ..., z[3]
  over rational field

julia> S_gens
Dict{Tuple{Union{Char, AbstractString, Symbol}, Edge}, MPolyRingElem} with 9 entries:
  (:x, Edge(4, 3)) => x[3]
  (:z, Edge(4, 3)) => z[3]
  (:z, Edge(4, 2)) => z[2]
  (:y, Edge(4, 1)) => y[1]
  (:x, Edge(4, 2)) => x[2]
  (:x, Edge(4, 1)) => x[1]
  (:y, Edge(4, 3)) => y[3]
  (:y, Edge(4, 2)) => y[2]
  (:z, Edge(4, 1)) => z[1]
```

### Model Ring and Equivalence Classes
The model ring's variables represent the joint probability distribution of states at the leaves of the tree. A key concept is that different leaf configurations can yield the same probability polynomial under the model's parametrization. These sets of configurations are called equivalence classes.

The full_model_ring contains a generator for every possible configuration of states at the leaves. A more compact and often more useful ring is the model_ring, which has one generator for each equivalence class.

```@docs
full_model_ring(PM::Union{PhylogeneticModel, GroupBasedPhylogeneticModel}; cached=false)
equivalent_classes(PM::Union{PhylogeneticModel, GroupBasedPhylogeneticModel})
model_ring(PM::Union{PhylogeneticModel, GroupBasedPhylogeneticModel}; cached=false)
```

## Parametrizations

A parametrization is a ring homomorphism from the model ring (representing observed probabilities) to the parameter ring (representing model parameters). This map explicitly gives the polynomial formula for each observable probability in terms of the model parameters.

There are two versions: `full_parametrization` maps from the `full_model_ring`, while `parametrization` maps from the reduced `model_ring`.

```@docs
full_parametrization(PM::PhylogeneticModel)
full_parametrization(PM::GroupBasedPhylogeneticModel)
parametrization(PM::PhylogeneticModel)
parametrization(PM::GroupBasedPhylogeneticModel)
```

### Affine Parametrization
Sometimes one wants to work on an affine variety by imposing constraints, such as the rows of the transition matrices summing to one. `affine_parametrization` provides these constrained maps.
```@docs
affine_parametrization(PM::PhylogeneticModel)
affine_parametrization(PM::GroupBasedPhylogeneticModel)
```

Example of a full workflow:

```jldoctest param
julia> tree = graph_from_edges(Directed,[[4,1], [4,2], [4,3]]);

julia> CFN = cavender_farris_neyman_model(tree);

# Get the reduced model ring and the parameter ring
julia> R, p = model_ring(phylogenetic_model(CFN));
julia> S, S_gens = parameter_ring(phylogenetic_model(CFN));

# Get the parametrization map
julia> f = parametrization(phylogenetic_model(CFN));

julia> f(p[(1, 1, 1)])  # Probability of observing state 1 at both leaves
1//2*a[1]*a[2]*b[3] + 1//2*a[3]*b[1]*b[2]
```


## Coordinate Change for Group-Based Models

For GroupBasedPhylogeneticModel objects, this module provides functions to convert between the standard probability coordinates and the Fourier coordinates. This change of basis often simplifies the phylogenetic invariants and the geometry of the model.

The functions fourier_transform and inverse_fourier_transform compute the transformation matrices, while coordinate_change and inverse_coordinate_change return the actual ring homomorphisms.

```@docs
fourier_transform(PM::GroupBasedPhylogeneticModel)
inverse_fourier_transform(PM::GroupBasedPhylogeneticModel)
coordinate_change(PM::GroupBasedPhylogeneticModel)
inverse_coordinate_change(PM::GroupBasedPhylogeneticModel)
```

```jldoctest change_coord
julia> tree = graph_from_edges(Directed,[[4,1], [4,2], [4,3]]);
julia> CFN = cavender_farris_neyman_model(tree);

# Get the coordinate change map from Fourier (q) to probability (p) coordinates
julia> cc = coordinate_change(CFN);

# Get the model rings for both coordinate systems
julia> Rq, q = model_ring(CFN);
julia> Rp, p = model_ring(phylogenetic_model(CFN));

# Express the first Fourier coordinate in terms of probabilities
julia> cc(q[(1, 1)])
p[1,2,2] + p[1,2,1] + p[1,1,2] + p[1,1,1]

julia> cc_inv(cc(q[(1, 1, 1)]))
q[1,1,1]
```


## Contact

Please direct questions about this part of OSCAR to the following people:
* [Marina Garrote López](https://sites.google.com/view/marinagarrotelopez)

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).
Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).
