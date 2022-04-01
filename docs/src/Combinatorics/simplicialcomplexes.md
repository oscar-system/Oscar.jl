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

### Subcomplexes

```@docs
star_subcomplex(K::SimplicialComplex, sigma::Union{Vector{Int}, Set{Int}})
link_subcomplex(K::SimplicialComplex, sigma::Union{Vector{Int}, Set{Int}})
```

### Surface examples

```@docs
torus()
klein_bottle()
real_projective_plane()
```

### Other examples

```@docs
complex_projective_plane()
```

## Basic properties

```@docs
nvertices(K::SimplicialComplex)
dim(K::SimplicialComplex)
f_vector(K::SimplicialComplex)
h_vector(K::SimplicialComplex)
euler_characteristic(K::SimplicialComplex)
```

## Homology and cohomology

```@docs
homology(K::SimplicialComplex, i::Int)
betti_numbers(K::SimplicialComplex)
cohomology(K::SimplicialComplex, i::Int)
```
## Fundamental group

```@docs
fundamental_group(K::SimplicialComplex)
```

## Connection to commutative algebra

The complements of the minimal non-faces form the facets of the Alexander dual.

```@docs
minimal_nonfaces(K::SimplicialComplex)
alexander_dual(K::SimplicialComplex)
```

Let $K$ be a simplicial complex on $n$ vertices.
The minimal non-faces of $K$ generate a square-free monomial ideal, known as the *Stanley-Reisner ideal* of $K$.
The quotient of the polynomial ring (in $n$ variables, with integer coefficients) modulo that ideal is the *Stanley-Reisner ring*.
For details see Chapter 5 of [BH09](@cite).

```@docs
stanley_reisner_ideal(K::SimplicialComplex)
stanley_reisner_ideal(R::MPolyRing, K::SimplicialComplex)
stanley_reisner_ring(K::SimplicialComplex)
stanley_reisner_ring(R::MPolyRing, K::SimplicialComplex)
```

## Saving and loading

Objects of type `SimplicialComplex` can be saved to a file and loaded with the
two methods `save` and `load`.  The file is in JSON format and contains the
underlying polymake object.  In particular, such a file can be read by both
polymake and Oscar.
