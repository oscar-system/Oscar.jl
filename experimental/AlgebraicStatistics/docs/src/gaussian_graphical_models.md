```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Gaussian Graphical Models

Gaussian graphical model types are parametrized as `GaussianGraphicalModel{T, L}` where `T` represents the type of graph and `L` represents the labelings of the graph. The parameters `T` and `L` are determined by the graph used to construct the model.

```@docs
gaussian_graphical_model
```

## Directed Gaussian Graphical Model 
The parametrization for a directed Gaussian graphical model on a DAG $G$ is built from the weighted adjacency matrix $\Lambda$ and the covariance matrix $\Omega$ of the error terms. The model is then the set of all covariance matrices $\Sigma = (Id - \Lambda)^{-T} \Omega (Id - \Lambda)^{-1}$. $\Lambda$ and $\Omega$ can be built with the following functions:
```@docs
directed_edges_matrix
error_covariance_matrix
```



## Undirected Gaussian Graphical Model
As with their directed counterpart, it is very easy to create new subtypes of graphical models by overloading the function `concentration_matrix` below. Unlike most other types of graphical models though, the vanishing ideal computation of an undirected graphical model is done by eliminating all concentration variables $k_{ij}$ from the ideal given by the equations $\Sigma K - Id$ after saturating by $\det(K)$. 


```@docs
concentration_matrix
```


## Gaussian Rings

```@docs
gaussian_ring
covariance_matrix
```
