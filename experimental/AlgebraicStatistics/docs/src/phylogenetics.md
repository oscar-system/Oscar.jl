```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Phylogenetics

The Phylogenetics module in OSCAR provides a comprehensive toolkit for the algebraic and geometric study of phylogenetic models, specifically nucleotide substitution Markov processes on phylogenetic trees. It allows for the construction of various evolutionary models on phylogenetic trees, the creation of their associated polynomial rings and parametrizations, and the computation of their algebraic properties, such as phylogenetic invariants (polynomials in the venishing ideal of the model).

## Model Constructors

Two main types of models can be constructed: a general `PhylogeneticModel` and a more specialized `GroupBasedPhylogeneticModel`.

#### `PhylogeneticModel`
A `PhylogeneticModel` is defined by a directed tree, a base field, the symbolic structure of the transition matrices for each edge, and a probability distribution at the root.

```@docs
PhylogeneticModel
phylogenetic_model
```

Here are some examples of how to construct phylogenetic models. The Jukes Cantor model on a 3-leaf tree can be constructed:
```jldoctest Example_F81:
julia> tree = graph_from_edges(Directed,[[4,1], [4,2], [4,3]]);

julia> M_JC = [:a :b :b :b;
               :b :a :b :b;
               :b :b :a :b;
               :b :b :b :a];

julia> PM_jukesCantor = phylogenetic_model(tree, M_JC)
Phylogenetic model on a tree with 3 leaves and 3 edges
with root distribution [1//4, 1//4, 1//4, 1//4] and transition matrices of the form
 [:a :b :b :b;
  :b :a :b :b;
  :b :b :a :b;
  :b :b :b :a].
```

And the Tamura-Nei model (TN93) in the same tree can be constructed as follows:
```jldoctest Example_F81:
julia> tree = graph_from_edges(Directed,[[4,1], [4,2], [4,3]]);

julia> Rp, r = rational_function_field(QQ, :r => 1:4);

julia> RM, (a, b, c, d) = polynomial_ring(Rp, [:a, :b, :c, :d]);

julia> M_TN93 = [a*r[1] c*r[2] b*r[3] b*r[4];
                c*r[1] a*r[2] b*r[3] b*r[4];
                b*r[1] b*r[2] a*r[3] d*r[4];
                b*r[1] b*r[2] d*r[3] a*r[4]];

julia> PM_TN93 = PhylogeneticModel(tree, M_TN93, r)
Phylogenetic model on a tree with 3 leaves and 3 edges
with root distribution [r[1], r[2], r[3], r[4]] and transition matrices of the form
 [r[1]*a r[2]*c r[3]*b r[4]*b;
  r[1]*c r[2]*a r[3]*b r[4]*b;
  r[1]*b r[2]*b r[3]*a r[4]*d;
  r[1]*b r[2]*b r[3]*d r[4]*a].
```

#### `GroupBasedPhylogeneticModel`

For models exhibiting symmetries that can be captured by a finite abelian group, the `GroupBasedPhylogeneticModel` is used. This structure requires the symbolic Fourier parameters and the associated group, in addition to the standard phylogenetic model components.

```@docs
GroupBasedPhylogeneticModel
group_based_phylogenetic_model
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
Group-based phylogenetic model on a tree with 3 leaves and 3 edges 
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
cavender_farris_neyman_model(F::Field, G::AbstractGraph{Directed})
jukes_cantor_model(F::Field, G::AbstractGraph{Directed})
kimura2_model(F::Field, G::AbstractGraph{Directed})
kimura3_model(F::Field, G::AbstractGraph{Directed})
general_markov_model(F::Field, G::AbstractGraph{Directed})
general_time_reversible_model(F::Field, G::AbstractGraph{Directed})
```

## Properties

Objects of type `GroupBasedPhylogeneticModel` are built upon a `PhylogeneticModel`, which can be accessed directly:
```@docs
phylogenetic_model(PM::GroupBasedPhylogeneticModel)
```

Other properties of a `PhylogeneticModel` and a `GroupBasedPhylogeneticModel` can be accessed using the following functions.

```@docs
base_field(PM::PhylogeneticModel)
n_states
root_distribution
transition_matrix
fourier_parameters
group(PM::GroupBasedPhylogeneticModel)
varnames(PM::PhylogeneticModel)
```

## Phylogenetic rings
Two main rings are associated with any phylogenetic model: the _parameter ring_, which contains the unobserved parameters of the evolutionary process (e.g., transition probabilities), and the _model ring_, which represents the observed data (e.g., joint probability distribution at the leaves).

### Parameter ring
The parameter ring's generators correspond to the symbolic parameters of the model. For a `PhylogeneticModel`, these are the entries of the transition matrices and root distribution. For a `GroupBasedPhylogeneticModel`, they are the Fourier parameters.

```@docs
parameter_ring(PM::PhylogeneticModel; cached=false, sorted_edges::Union{Vector{Edge}, Nothing} = nothing)
parameter_ring(PM::GroupBasedPhylogeneticModel; cached=false, sorted_edges::Union{Vector{Edge}, Nothing} = nothing)
```

To access specific parameters (entries of the treansition matrices, fourier parameters, entry of the root distriution or paramters associated to hybrid nodes in phylo networks) you can use the following functions:

```@docs
entry_root_distribution
entry_transition_matrix
entry_fourier_parameter
entry_hybrid_parameter
```

For example

```jldoctest Example_JC:
julia> tree = graph_from_edges(Directed,[[4,1], [4,2], [4,3]]);

julia> PM = kimura2_model(tree)
Group-based phylogenetic model on a tree with 3 leaves and 3 edges
with root distribution [1//4, 1//4, 1//4, 1//4],
transition matrices of the form
 [:a :b :c :b;
  :b :a :b :c;
  :c :b :a :b;
  :b :c :b :a]
and fourier parameters of the form [:x, :y, :z, :z].

julia> S, - = parameter_ring(PM);

julia> S
Multivariate polynomial ring in 9 variables x[1], x[2], x[3], y[1], ..., z[3]
  over rational field

julia> entry_root_distribution(PM, 1)
1//4

julia> entry_transition_matrix(PM, 3, 3, Edge(4, 1))
a[1]

julia> entry_fourier_parameter(PM, 4, Edge(4, 1))
z[1]
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
full_parametrization
parametrization(PM::Union{PhylogeneticModel, GroupBasedPhylogeneticModel})
```

*Affine Parametrization:* 
Sometimes one wants to work on an affine variety by imposing constraints, such as the rows of the transition matrices summing to one. `affine_parametrization` provides these constrained maps.
```@docs
affine_parametrization
```

Example of a full workflow:

```jldoctest param
julia> tree = graph_from_edges(Directed,[[4,1], [4,2], [4,3]]);

julia> JC = jukes_cantor_model(tree);

julia> pmJC = phylogenetic_model(JC);

julia> Rp, p = model_ring(pmJC);

julia> Rq, q = model_ring(JC);

julia> f = parametrization(pmJC);

julia> f(p[(1, 1, 1)])  # Probability of observing state 1 at all leaves
1//4*a[1]*a[2]*a[3] + 3//4*b[1]*b[2]*b[3]

julia> g = parametrization(JC);

julia> g(q[(1, 1, 1)]) 
x[1]*x[2]*x[3]

julia> vanishing_ideal(pmJC)
Ideal generated by
  2*p[1,2,3]^3 - p[1,2,3]^2*p[1,2,2] - p[1,2,3]^2*p[1,2,1] - p[1,2,3]^2*p[1,1,2] + p[1,2,3]^2*p[1,1,1] - p[1,2,3]*p[1,2,2]^2 - p[1,2,3]*p[1,2,1]^2 - p[1,2,3]*p[1,1,2]^2 + p[1,2,3]*p[1,1,1]^2 + p[1,2,2]^2*p[1,2,1] + p[1,2,2]^2*p[1,1,2] + p[1,2,2]*p[1,2,1]^2 - p[1,2,2]*p[1,2,1]*p[1,1,2] - p[1,2,2]*p[1,2,1]*p[1,1,1] + p[1,2,2]*p[1,1,2]^2 - p[1,2,2]*p[1,1,2]*p[1,1,1] + p[1,2,1]^2*p[1,1,2] + p[1,2,1]*p[1,1,2]^2 - p[1,2,1]*p[1,1,2]*p[1,1,1]

julia> vanishing_ideal(JC)
Ideal generated by
  -q[2,3,4]^2*q[1,1,1] + q[2,2,1]*q[2,1,2]*q[1,2,2]
```

The vanishing ideal can be computed using different algorithms

```jldoctest param
julia> tree = graph_from_edges(Directed,[[4,1], [4,2], [4,3]]);

julia> JC = jukes_cantor_model(tree);

julia> vanishing_ideal(JC; algorithm = :kernel)
Ideal generated by
  -q[2,3,4]^2*q[1,1,1] + q[2,2,1]*q[2,1,2]*q[1,2,2]

julia> vanishing_ideal(JC; algorithm = :eliminate)
Ideal generated by
  -q[2,3,4]^2*q[1,1,1] + q[2,2,1]*q[2,1,2]*q[1,2,2]

julia> vanishing_ideal(JC; algorithm = :f4)
Ideal generated by
  -q[2,3,4]^2*q[1,1,1] + q[2,2,1]*q[2,1,2]*q[1,2,2]

julia> vanishing_ideal(JC; algorithm = :markov)
Ideal generated by
  q[2,3,4]^2*q[1,1,1] - q[2,2,1]*q[2,1,2]*q[1,2,2]
```

## Coordinate Change for Group-Based Models

For `GroupBasedPhylogeneticModel` objects, this module provides functions to convert between the standard probability coordinates and the Fourier coordinates. This change of basis often simplifies the phylogenetic invariants and the geometry of the model.

The functions `fourier_transform` and `inverse_fourier_transform` compute the transformation matrices, while `coordinate_change` and `inverse_coordinate_change` return the actual ring homomorphisms.

```@docs
fourier_transform(PM::GroupBasedPhylogeneticModel)
inverse_fourier_transform(PM::GroupBasedPhylogeneticModel)
coordinate_change(PM::GroupBasedPhylogeneticModel)
inverse_coordinate_change(PM::GroupBasedPhylogeneticModel)
```

```jldoctest change_coord
julia> tree = graph_from_edges(Directed,[[4,1], [4,2], [4,3]]);

julia> CFN = cavender_farris_neyman_model(tree);

julia> cc = coordinate_change(CFN);

julia> cc_inv = inverse_coordinate_change(CFN);

julia> Rq, q = model_ring(CFN);

julia> Rp, p = model_ring(phylogenetic_model(CFN));

julia> cc(q[(1, 1, 1)])
p[1,2,2] + p[1,2,1] + p[1,1,2] + p[1,1,1]

julia> cc_inv(cc(q[(1, 1, 1)]))
q[1,1,1]
```


## Phylogenetic Networks

The `PhylogeneticNetwork` structure provides a way to represent and manipulate **phylogenetic networks**, which are directed acyclic graphs (DAGs) that generalizes phylogenetic trees by allowing hybrid nodes — nodes with multiple incoming edges — that represent genetic mixing events between lineages.

```@docs
PhylogeneticNetwork
phylogenetic_network
```

You can access specific properties of a `PhylogeneticNetwork`, like its level and hybrid nodes, using these functions:

```@docs
level(N::PhylogeneticNetwork)
n_hybrid
hybrids
hybrid_vertices
hybrid_edges
```

General graph properties can be accessed with the following functions:
```@docs
graph(N::PhylogeneticNetwork)
n_edges(N::PhylogeneticNetwork)
edges(N::PhylogeneticNetwork)
n_leaves(N::PhylogeneticNetwork)
leaves(N::PhylogeneticNetwork)
```

All model constructors, including `PhylogeneticModel`, `GroupBasedPhylogeneticModel`, and the classic pre-defined models, can be applied to phylogenetic networks in the same way they are applied to trees. The resulting model object will be defined on the network, incorporating its specific topology. Consequently, all network property functions (e.g., graph, level, hybrids) can be called directly on the model object.

### Example: CFN model on a Level-1 Network
The following example demonstrates this by creating a Cavender-Farris-Neyman model on a level-1 network. It shows how to inspect the model's underlying graph and access the algebraic parameters. 
We use the method `components_of_kernel(d::Int, phi::MPolyAnyMap)` with $d=3$ to compute all minimal generators of the kernel of the polynomial map $\phi$ with total degree at most $3$ using the main algorithm of [CH26](@cite).


```jldoctest phylo_network
julia> G = graph_from_edges(Directed,[[7,6], [7,8], [6,5], [8,5], [5,1], [6,2], [7,3], [8,4]]);

julia> PM = cavender_farris_neyman_model(G)
Group-based phylogenetic model on a level-1 network with 1 hybrid node, 4 leaves  
and 8 edges with root distribution [1//2, 1//2], 
transition matrices of the form 
 [:a :b;
  :b :a]
and fourier parameters of the form [:x, :y].

julia> graph(PM)
Level-1 phylogenetic network with hybrid nodes {5} and edges
  (5, 1)(6, 2)(6, 5)(7, 3)(7, 6)(7, 8)(8, 4)(8, 5)

julia> S, - = parameter_ring(PM);

julia> S
Multivariate polynomial ring in 18 variables l[1, 1], l[1, 2], x[1], x[2], ..., y[8]
  over rational field

julia> entry_fourier_parameter(PM, 1, Edge(6, 5))
x[5]

julia> entry_hybrid_parameter(PM, Edge(6, 5))
l[1, 1]

julia> R, q = full_model_ring(PM);

julia> f = full_parametrization(PM);

julia> f(q[1,1,1,1])
l[1, 1]*x[1]*x[2]*x[3]*x[4]*x[5]*x[6]*x[7] + l[1, 2]*x[1]*x[2]*x[3]*x[4]*x[6]*x[7]*x[8]

julia> phi = parametrization(PM)
Ring homomorphism
  from multivariate polynomial ring in 8 variables over QQ
  to multivariate polynomial ring in 18 variables over QQ
defined by
  q[2,2,2,2] -> l[1, 1]*x[6]*y[1]*y[2]*y[3]*y[4]*y[5]*y[7] + l[1, 2]*x[7]*y[1]*y[2]*y[3]*y[4]*y[6]*y[8]
  q[2,2,1,1] -> l[1, 1]*x[3]*x[4]*x[6]*x[7]*y[1]*y[2]*y[5] + l[1, 2]*x[3]*x[4]*y[1]*y[2]*y[6]*y[7]*y[8]
  q[2,1,2,1] -> l[1, 1]*x[2]*x[4]*x[7]*y[1]*y[3]*y[5]*y[6] + l[1, 2]*x[2]*x[4]*x[6]*y[1]*y[3]*y[7]*y[8]
  q[2,1,1,2] -> l[1, 1]*x[2]*x[3]*y[1]*y[4]*y[5]*y[6]*y[7] + l[1, 2]*x[2]*x[3]*x[6]*x[7]*y[1]*y[4]*y[8]
  q[1,2,2,1] -> l[1, 1]*x[1]*x[4]*x[5]*x[7]*y[2]*y[3]*y[6] + l[1, 2]*x[1]*x[4]*x[7]*x[8]*y[2]*y[3]*y[6]
  q[1,2,1,2] -> l[1, 1]*x[1]*x[3]*x[5]*y[2]*y[4]*y[6]*y[7] + l[1, 2]*x[1]*x[3]*x[8]*y[2]*y[4]*y[6]*y[7]
  q[1,1,2,2] -> l[1, 1]*x[1]*x[2]*x[5]*x[6]*y[3]*y[4]*y[7] + l[1, 2]*x[1]*x[2]*x[6]*x[8]*y[3]*y[4]*y[7]
  q[1,1,1,1] -> l[1, 1]*x[1]*x[2]*x[3]*x[4]*x[5]*x[6]*x[7] + l[1, 2]*x[1]*x[2]*x[3]*x[4]*x[6]*x[7]*x[8]

julia> collect(values(components_of_kernel(3, phi)))
9-element Vector{Vector{QQMPolyRingElem}}:
 [q[2,2,2,2]*q[2,1,2,1]*q[1,1,1,1] - q[2,2,1,1]*q[2,1,2,1]*q[1,1,2,2] + q[2,1,2,1]^2*q[1,2,1,2] - q[2,1,2,1]*q[2,1,1,2]*q[1,2,2,1]]
 [q[2,2,2,2]*q[1,1,1,1] - q[2,2,1,1]*q[1,1,2,2] + q[2,1,2,1]*q[1,2,1,2] - q[2,1,1,2]*q[1,2,2,1]]
 [q[2,2,2,2]*q[1,1,1,1]^2 - q[2,2,1,1]*q[1,1,2,2]*q[1,1,1,1] + q[2,1,2,1]*q[1,2,1,2]*q[1,1,1,1] - q[2,1,1,2]*q[1,2,2,1]*q[1,1,1,1]]
 [q[2,2,2,2]*q[1,2,2,1]*q[1,1,1,1] - q[2,2,1,1]*q[1,2,2,1]*q[1,1,2,2] + q[2,1,2,1]*q[1,2,2,1]*q[1,2,1,2] - q[2,1,1,2]*q[1,2,2,1]^2]
 [q[2,2,2,2]*q[2,1,1,2]*q[1,1,1,1] - q[2,2,1,1]*q[2,1,1,2]*q[1,1,2,2] + q[2,1,2,1]*q[2,1,1,2]*q[1,2,1,2] - q[2,1,1,2]^2*q[1,2,2,1]]
 [q[2,2,2,2]*q[2,2,1,1]*q[1,1,1,1] - q[2,2,1,1]^2*q[1,1,2,2] + q[2,2,1,1]*q[2,1,2,1]*q[1,2,1,2] - q[2,2,1,1]*q[2,1,1,2]*q[1,2,2,1]]
 [q[2,2,2,2]*q[1,1,2,2]*q[1,1,1,1] - q[2,2,1,1]*q[1,1,2,2]^2 + q[2,1,2,1]*q[1,2,1,2]*q[1,1,2,2] - q[2,1,1,2]*q[1,2,2,1]*q[1,1,2,2]]
 [q[2,2,2,2]^2*q[1,1,1,1] - q[2,2,2,2]*q[2,2,1,1]*q[1,1,2,2] + q[2,2,2,2]*q[2,1,2,1]*q[1,2,1,2] - q[2,2,2,2]*q[2,1,1,2]*q[1,2,2,1]]
 [q[2,2,2,2]*q[1,2,1,2]*q[1,1,1,1] - q[2,2,1,1]*q[1,2,1,2]*q[1,1,2,2] + q[2,1,2,1]*q[1,2,1,2]^2 - q[2,1,1,2]*q[1,2,2,1]*q[1,2,1,2]]

julia> coordinate_change(PM)
Ring homomorphism
  from multivariate polynomial ring in 8 variables over QQ
  to multivariate polynomial ring in 8 variables over QQ
defined by
  q[2,2,2,2] -> -p[1,2,2,2] + p[1,2,2,1] + p[1,2,1,2] - p[1,2,1,1] + p[1,1,2,2] - p[1,1,2,1] - p[1,1,1,2] + p[1,1,1,1]
  q[2,2,1,1] -> -p[1,2,2,2] - p[1,2,2,1] - p[1,2,1,2] - p[1,2,1,1] + p[1,1,2,2] + p[1,1,2,1] + p[1,1,1,2] + p[1,1,1,1]
  q[2,1,2,1] -> -p[1,2,2,2] - p[1,2,2,1] + p[1,2,1,2] + p[1,2,1,1] - p[1,1,2,2] - p[1,1,2,1] + p[1,1,1,2] + p[1,1,1,1]
  q[2,1,1,2] -> -p[1,2,2,2] + p[1,2,2,1] - p[1,2,1,2] + p[1,2,1,1] - p[1,1,2,2] + p[1,1,2,1] - p[1,1,1,2] + p[1,1,1,1]
  q[1,2,2,1] -> p[1,2,2,2] + p[1,2,2,1] - p[1,2,1,2] - p[1,2,1,1] - p[1,1,2,2] - p[1,1,2,1] + p[1,1,1,2] + p[1,1,1,1]
  q[1,2,1,2] -> p[1,2,2,2] - p[1,2,2,1] + p[1,2,1,2] - p[1,2,1,1] - p[1,1,2,2] + p[1,1,2,1] - p[1,1,1,2] + p[1,1,1,1]
  q[1,1,2,2] -> p[1,2,2,2] - p[1,2,2,1] - p[1,2,1,2] + p[1,2,1,1] + p[1,1,2,2] - p[1,1,2,1] - p[1,1,1,2] + p[1,1,1,1]
  q[1,1,1,1] -> p[1,2,2,2] + p[1,2,2,1] + p[1,2,1,2] + p[1,2,1,1] + p[1,1,2,2] + p[1,1,2,1] + p[1,1,1,2] + p[1,1,1,1]


```
