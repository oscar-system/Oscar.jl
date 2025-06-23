# Gaussian Graphical Models

The OSCAR type for graphical models is of parametrized form `GraphicalModel{G, T}` where `T` represents the type of ring in which the vanishing ideal of the model belongs and `G` represents the associated graph. Gaussian graphical models are those where the `T` is of type `GaussianRing` which is a multivariate polynomial ring equipped with some extra features.

## Gaussian Rings

```@docs
gaussian_ring(n::Int; s_var_name::VarName="s", K::Field=QQ, cached=false)
gens(R::GaussianRing)
covariance_matrix(R::GaussianRing)
```

## Directed Gaussian Graphical Models

A directed Gaussian graphical model is constructed from `G::Graph{Directed}` and `S::GaussianRing`. Optionally, the user may specify a string `l_var_name::String` which corresponds to the edge weights in the parametrization of the model and a string `w_var_name::String` for labeling the error covariance parameters.
```@docs
graphical_model(G::Graph{Directed}, S::GaussianRing; l_var_name::VarName="l", w_var_name::VarName="w", cached=false)
```

The parametrization for a directed Gaussian graphical model on a DAG $G$ is built from the weighted adjacency matrix $\Lambda$ and the covariance matrix $\Omega$ of the error terms. The model is then the set of all covariance matrices $\Sigma = (Id - \Lambda)^{-T} \Omega (Id - \Lambda)^{-1}$. $\Lambda$ and $\Omega$ can be built with the following functions:
```@docs
directed_edges_matrix(M::GraphicalModel{Graph{Directed}, GaussianRing})
error_covariance_matrix(M::GraphicalModel{Graph{Directed}, GaussianRing})
```

It is easy to create new types of graphical models by overloading the methods `directed_edges_matrix` and `error_covariance_matrix`. For instance, colored graphical models may easily be created by creating a new type of `GraphicalModel{G, T}` and a new `directed_edges_matrix` function. The following two functions will then work almost immediately on this new type.
```@docs
parametrization(M::GraphicalModel{Graph{Directed}, GaussianRing})
vanishing_ideal(M::GraphicalModel{Graph{Directed}, GaussianRing})
```

With almost all graphical models, the vanishing ideal is computed by taking the kernel of the ring map given by `parametrization(M:GraphicalModel{G, T})`. 



## Undirected Gaussian Graphical Models


A undirected Gaussian graphical model is constructed from `G::Graph{Undirected}`, `S::GaussianRing`. Optionally, the user may specify a string `k_var_name::String` which corresponds to the entries of the concentration matrix.
```@docs
graphical_model(G::Graph{Undirected}, S::GaussianRing; k_var_name::VarName="k", cached=false)
```

As with their directed counterpart, it is very easy to create new subtypes of graphical models by overloading the function `concentration_matrix` below. Unlike most other types of graphical models though, the vanishing ideal computation of an undirected graphical model is done by eliminating all concentration variables $k_{ij}$ from the ideal given by the equations $\Sigma K - Id$ after saturating by $\det(K)$. 


```@docs
concentration_matrix(M::GraphicalModel{Graph{Undirected}, GaussianRing})
parametrization(M::GraphicalModel{Graph{Undirected}, GaussianRing})
vanishing_ideal(M::GraphicalModel{Graph{Undirected}, GaussianRing})
```
