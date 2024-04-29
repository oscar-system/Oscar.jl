```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```
# Algebraic Phylogenetics

Add some maths introduction here.

The Algebraic Phylogenetics part of OSCAR provides functionality for
- specifying phylogenetic models
- parametrizing such models
- calculating their algebraic invariants

## Entities and their components

This package provides functionality for two basic structures: `PhylogeneticModel` specifying a phylogenetic tree model in probability coordinates and `GroupBasedPhylogeneticModel` for group based, e.g. Fourier, coordinates. Such a model is made up of a graph, transition matrices attached to the graph's edges, the number of states of each vertex random variable, and the corresponding polynomial rings.

```@docs
  graph(pm::PhylogeneticModel)
  transition_matrices(pm::PhylogeneticModel)
  number_states(pm::PhylogeneticModel)
  probability_ring(pm::PhylogeneticModel)
  fourier_ring(pm::PhylogeneticModel)
  fourier_parameters(pm::PhylogeneticModel)
  group_model(pm::PhylogeneticModel)
```

## Models

The five most common models in algebraic phylogenetics can automatically be specified by calling the below functions, each taking a tree as input and attaching transition matrices to its edges as defined by Jukes Cantor, Kimura, etc. respectively, returning a structure `PhylogeneticModel` or `GroupBasedPhylogeneticModel`.

```@docs
  jukes_cantor_model(graph::Graph{Directed})
  kimura2_model(graph::Graph{Directed})
  kimura3_model(graph::Graph{Directed})
  general_markov_model(graph::Graph{Directed})
  cavender_farris_neyman_model(graph::Graph{Directed})
```

## Parametrisation

For each phylogenetic model, we can calculate the parametrisation, a map from transition matrices to probabilities, parametrized in probability or Fourier coordinates. We can reparametrize between these and return the transformation matrix, and we can calculate equivalence classes of probabilities with the same parametrization.

```@docs
  probability_map(pm::PhylogeneticModel)
  fourier_map(pm::GroupBasedPhylogeneticModel)
  specialized_fourier_transform(pm::GroupBasedPhylogeneticModel, p_equivclasses::Dict{Vector{Vector{Int64}}, QQMPolyRingElem}, f_equivclasses::Dict{Vector{Vector{Int64}}, QQMPolyRingElem})
  inverse_specialized_fourier_transform(pm::GroupBasedPhylogeneticModel, p_equivclasses::Dict{Vector{Vector{Int64}}, QQMPolyRingElem},f_equivclasses::Dict{Vector{Vector{Int64}}, QQMPolyRingElem})

  compute_equivalent_classes
  sum_equivalent_classes(pm::PhylogeneticModel, equivalent_classes::Dict{Vector{Vector{Int64}}, QQMPolyRingElem})
```

## Contact

Please direct questions about this part of OSCAR to the following people:
* [Marina Garrote Lopéz](s://sites.google.com/view/marinagarrotelopez),
* possibly others.

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).
Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).
