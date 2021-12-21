```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["simplicialcomplexes.md"]
```

# Simplicial Complexes

## Introduction

Abstract simplicial complexes provide a combinatorial way to define topological spaces.
By no means every topological space arises in this way, but this is a (most) natural choice in a computational setup.

A *simplicial complex* $K$ on a vertex set $V$ is a nonempty subset of $2^V$ such that for each $\sigma \in K$ and $\tau \subset\sigma$ we have $\tau\in K$.
Here $V$ is usually $[n] = \{1,2,\dots,n\}$ for some $n\geq 0$.

General textbooks offering details on the theory include:
- [Koz08](@cite)


## Construction

```@docs
SimplicialComplex(K::Vector{Vector{Int}})
```

### Examples
```@docs
torus()
klein_bottle()
realprojectiveplane()
complexprojectiveplane()
```

## Basic properties
```@docs
nvertices(K::SimplicialComplex)
dim(K::SimplicialComplex)
f_vector(K::SimplicialComplex)
h_vector(K::SimplicialComplex)
```

## Connection to commutative algebra
```@docs
minimalnonfaces(K::SimplicialComplex)
stanley_reisner_ideal(K::SimplicialComplex)
stanley_reisner_ring(K::SimplicialComplex)
```

## Saving and loading

Objects of type `SimplicialComplex` can be saved to a file and loaded with the following two methods:
```@docs
save_simplicialcomplex(K::SimplicialComplex, filename::String)
load_simplicialcomplex(filename::String)
```
The file is in JSON format and contains the underlying polymake object.
In particular, such a file can be read by both polymake and Oscar.

