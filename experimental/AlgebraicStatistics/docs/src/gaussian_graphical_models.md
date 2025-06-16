```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Gaussian Graphical Models

Gaussian graphical model types are parametrized as `GaussianGraphicalModel{T, L} <: GraphicalModel{T, L}` where `T` represents the type of graph and `L` represents the labelings of the graph. 

```@docs
gaussian_graphical_model
```

The parametrization for a directed Gaussian graphical model on a DAG $G$ is built from the weighted adjacency matrix $\Lambda$ and the covariance matrix $\Omega$ of the error terms. The model is then the set of all covariance matrices $\Sigma = (Id - \Lambda)^{-T} \Omega (Id - \Lambda)^{-1}$. $\Lambda$ and $\Omega$ can be built with the following functions:
```@docs
directed_edges_matrix(M::GaussianGraphicalModel{Directed, L}) where L
error_covariance_matrix(M::GaussianGraphicalModel{Directed, L}) where L
```

With almost all graphical models, the vanishing ideal is computed by taking the kernel of the ring map given by `parametrization(M::GaussianGraphicalModel)`. 


## Undirected Gaussian Graphical Models


A undirected Gaussian graphical model is constructed from `G::Graph{Undirected}`. Optionally, the user may specify a string `k_var_name::String` which corresponds to the entries of the concentration matrix.
```@docs
gaussian_graphical_model(G::Graph{Undirected}; k_var_name::VarName="k", cached=false)
```

As with their directed counterpart, it is very easy to create new subtypes of graphical models by overloading the function `concentration_matrix` below. Unlike most other types of graphical models though, the vanishing ideal computation of an undirected graphical model is done by eliminating all concentration variables $k_{ij}$ from the ideal given by the equations $\Sigma K - Id$ after saturating by $\det(K)$. 


```@docs
concentration_matrix(M::GaussianGraphicalModel{Undirected, L}) where L
parametrization(M::GraphicalModel{Undirected, L}) where L
vanishing_ideal(M::GraphicalModel{Undirected, L}) where L
```

## Gaussian Rings

```@docs
gaussian_ring
covariance_matrix
```
