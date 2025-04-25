```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Algebraic Phylogenetics

The Algebraic Phylogenetics part of OSCAR provides functionality for
- specifying phylogenetic models
- parametrizing such models
- calculating their algebraic invariants

## Models

The five most common models in algebraic phylogenetics can automatically be specified by calling the below functions, each taking a tree as input and attaching transition matrices to its edges as defined by Jukes Cantor, Kimura, etc., respectively, returning a structure `PhylogeneticModel` or `GroupBasedPhylogeneticModel`.

```@docs
cavender_farris_neyman_model(graph::Graph{Directed})
jukes_cantor_model(graph::Graph{Directed})
kimura2_model(graph::Graph{Directed})
kimura3_model(graph::Graph{Directed})
general_markov_model(graph::Graph{Directed})
```

The models are by default in projective space. For their affine versions, call

```@docs
affine_phylogenetic_model!(pm::PhylogeneticModel)
```

## Components of a model

`PhylogeneticModel` specifies any phylogenetic tree model in probability coordinates and `GroupBasedPhylogeneticModel` can specify group-based, e.g. Fourier, coordinates. For any model, we can call its graph, transition matrices attached to the graph's edges, the number of states of each vertex random variable, and the corresponding polynomial rings. For instance for Jukes Cantor on the star with three leaves:

```@docs
graph(pm::PhylogeneticModel)
transition_matrices(pm::PhylogeneticModel)
number_states(pm::PhylogeneticModel)
probability_ring(pm::PhylogeneticModel)
fourier_ring(pm::GroupBasedPhylogeneticModel)
fourier_parameters(pm::GroupBasedPhylogeneticModel)
group_of_model(pm::GroupBasedPhylogeneticModel)
```


## Parametrization

For each phylogenetic model, we can calculate the parametrization, a map from transition matrices to probabilities, parametrized in probability or Fourier coordinates. For group-based models, we can reparametrize between these and return the transformation matrix, and we can calculate equivalence classes of probabilities with the same parametrization.

```@docs
probability_map(pm::PhylogeneticModel)
fourier_map(pm::GroupBasedPhylogeneticModel)
compute_equivalent_classes(parametrization::Dict{Tuple{Vararg{Int64}}, QQMPolyRingElem})
sum_equivalent_classes(equivalent_classes::NamedTuple{(:parametrization, :classes), Tuple{Dict{Tuple{Vararg{Int64}}, QQMPolyRingElem}, Dict{Tuple{Vararg{Int64}}, Vector{Tuple{Vararg{Int64}}}}}})
specialized_fourier_transform(pm::GroupBasedPhylogeneticModel, p_classes::Dict{Tuple{Vararg{Int64}}, Vector{Tuple{Vararg{Int64}}}}, q_classes::Dict{Tuple{Vararg{Int64}}, Vector{Tuple{Vararg{Int64}}}})
inverse_specialized_fourier_transform(pm::GroupBasedPhylogeneticModel, p_classes::Dict{Tuple{Vararg{Int64}}, Vector{Tuple{Vararg{Int64}}}}, q_classes::Dict{Tuple{Vararg{Int64}}, Vector{Tuple{Vararg{Int64}}}})
```

## Contact

Please direct questions about this part of OSCAR to the following people:
* [Marina Garrote López](https://sites.google.com/view/marinagarrotelopez),
* possibly others.

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).
Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).
