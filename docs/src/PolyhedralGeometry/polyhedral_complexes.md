```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["polyhedral_complexes.md"]
```

# Polyhedral Complexes

## Introduction

Let $\mathbb{F}$ be an ordered field; the default is that
$\mathbb{F}=\mathbb{Q}$ is the field of rational numbers and other fields are
not yet supported everywhere in the implementation.

A nonempty finite collection $\mathcal{P}$ of polyhedra in
$\mathbb{F}^n$, for $n$ fixed, is a *polyhedral complex* if

- the set $\mathcal{F}$ is closed with respect to taking faces and
- if $C,D\in\mathcal{F}$ then $C\cap D$ is a face of both, $C$ and $D$.

## Construction

To construct a polyhedral complex, you must pass points of each polyhedron in
the polyhedral complex, such that the polyhedron is the convex hull thereof,
along with an `IncidenceMatrix` encoding which points generate which
polyhedron.

```@docs
PolyhedralComplex
```


## Auxiliary functions
```@docs
ambient_dim(PC::PolyhedralComplex)
maximal_cones(PC::PolyhedralComplex)
nmaximal_cones(PC::PolyhedralComplex)
vertices(PC::PolyhedralComplex)
```

