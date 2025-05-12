```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Graphical Models

The OSCAR type for graphical models is of parametrized form `GraphicalModel{G, T}` where `T` represents the type of ring in which the vanishing ideal of the model belongs and `G` represents the associated graph. This parametrized typing allows the user to easily build upon the existing functionality to work with newer variations on graphical models such as those with colours or hidden variables. 

Many basic functions are defined for all graphical models `M` such as

- `graph(M)` refers to the associated graph
- `ring(M)` to the multivariate polynomial ring where the model resides
- `param_ring` to the multivariate polynomial ring where the parameters reside
- `param_gens` to the parameters of the model which are ring variables stored typically stored in a hash table for convenient indexing

This part of OSCAR also includes some basic graph functionality such as finding the vertices of a graph or its set of maximal cliques. Lastly, while many methods for graphical models depend heavily on whether or not they are discrete/Gaussian or directed/undirected some functionality is independent of this and thus implemented simultaneously for all graphical models. For example, the vanishing ideal of the model:

```@docs
vanishing_ideal(M::GraphicalModel)
```
