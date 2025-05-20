```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Polyhedral Complexes

## Introduction

Let $\mathbb{F}$ be an ordered field; the default is that
$\mathbb{F}=\mathbb{Q}$ is the field of rational numbers and other fields are
not yet supported everywhere in the implementation.

A nonempty finite collection $\mathcal{P}$ of polyhedra in
$\mathbb{F}^n$, for $n$ fixed, is a *polyhedral complex* if

- the set $\mathcal{F}$ is closed with respect to taking faces and
- if $C,D\in\mathcal{F}$ then $C\cap D$ is a face of both $C$ and $D$.

## Construction

To construct a polyhedral complex, you must pass points of each polyhedron in
the polyhedral complex, such that the polyhedron is the convex hull thereof,
along with an `IncidenceMatrix` encoding which points generate which
polyhedron.

```@docs
polyhedral_complex
```


## Auxiliary functions
```@docs
ambient_dim(PC::PolyhedralComplex)
codim(PC::PolyhedralComplex)
dim(PC::PolyhedralComplex)
f_vector(PC::PolyhedralComplex)
is_embedded(PC::PolyhedralComplex)
is_pure(PC::PolyhedralComplex)
is_simplicial(PC::PolyhedralComplex)
lineality_dim(PC::PolyhedralComplex)
lineality_space(PC::PolyhedralComplex{T}) where T<:scalar_types
maximal_polyhedra(PC::PolyhedralComplex{T}) where T<:scalar_types
minimal_faces(PC::PolyhedralComplex{T}) where T<:scalar_types
n_maximal_polyhedra(PC::PolyhedralComplex)
n_polyhedra(PC::PolyhedralComplex)
n_rays(PC::PolyhedralComplex)
n_vertices(PC::PolyhedralComplex)
polyhedra_of_dim
rays(PC::PolyhedralComplex{T}) where T<:scalar_types
rays_modulo_lineality(PC::PolyhedralComplex{T}) where T<:scalar_types
vertices(as::Type{PointVector{T}}, PC::PolyhedralComplex{T}) where {T<:scalar_types}
vertices_and_rays(PC::PolyhedralComplex{T}) where T<:scalar_types
```
