```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["cones.md"]
```

# Cones


## Introduction

Let $\mathbb{F}$ be an ordered field; the default is that
$\mathbb{F}=\mathbb{Q}$ is the field of rational numbers and other fields are
not yet supported everywhere in the implementation.

A set $C \subseteq \mathbb{F}^n$ is called a *(polyhedral) cone* if it can be
written as the set of nonnegative linear combinations of finitely many vectors
in $\mathbb{F}^n$.  Equivalently, cones can be written as the intersection of
finitely many homogeneous linear inequalities.

Any cone is a special case of a polyhedron.  Conversely, intersecting a cone
with a suitable affine hyperplane yields a polyhedron whose faces are in
bijection with the faces of the cone.  Going back and forth between polyhedra
and their homogenizations, the cones, is a frequent operation.  This is one
reason for keeping cones as a distinct type.

## Construction

```@docs
positive_hull(::Type{T}, ::Union{Oscar.MatElem, AbstractMatrix}) where T<:scalar_types
secondary_cone(SOP::SubdivisionOfPoints{T}) where T<:scalar_types
```

## Saving and loading

Objects of type `Cone` can be saved to a file and loaded from a file in the
following way:
```@repl oscar
C = positive_hull([1 0; 0 1])
save(C, "C.cone")
CC = load("C.cone")
collect(rays(CC))
```
The file is in JSON format and contains all previously gathered data belonging
to the underlying polymake object. In particular, this file can now be read by
both polymake and Oscar.

## Auxiliary functions
```@docs
ambient_dim(C::Cone)
contains(C::Cone, v::AbstractVector)
f_vector(C::Cone)
hilbert_basis(C::Cone{fmpq})
codim(C::Cone)
dim(C::Cone)
polarize(C::Cone{T}) where T<:scalar_types
intersect(C0::Cone{T}, C1::Cone{T}) where T<:scalar_types
is_pointed(C::Cone)
is_fulldimensional(C::Cone)
lineality_dim(C::Cone)
lineality_space(C::Cone{T}) where T<:scalar_types
nfacets(C::Cone)
nrays(C::Cone)
rays(C::Cone)
```

### Visualization
```@docs
visualize(C::Cone)
```
