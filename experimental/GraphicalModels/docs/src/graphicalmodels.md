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

`PhylogeneticModel` for probability coordinates, `FourierPhylogeneticModel` for Fourier coordinates.

```@docs
  graph(pm::PhylogeneticModel)
  transition_matrices(pm::PhylogeneticModel)
  number_states(pm::PhylogeneticModel)
  probabilities_ring(pm::PhylogeneticModel)

  fourier_ring(pm::PhylogeneticModel)
  fourier_parameters(pm::PhylogeneticModel)
  group_model(pm::PhylogeneticModel)
```

## Models

```@docs
  jukes_cantor_model(graph::Graph{Directed})
  kimura2_model(graph::Graph{Directed})
  kimura3_model(graph::Graph{Directed})
  generalmarkov_model(graph::Graph{Directed})
```


## Contact

Please direct questions about this part of OSCAR to the following people:
* [Marina Garrote Lopéz](s://sites.google.com/view/marinagarrotelopez),
* possibly others.

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).
Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).
