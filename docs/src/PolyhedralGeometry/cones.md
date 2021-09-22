```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["pg_cones.md"]
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
positive_hull(::Union{Oscar.MatElem, AbstractMatrix})
secondary_cone(SOP::SubdivisionOfPoints)
```

## Saving and loading

Objects of type `Cone` can be saved to a file and loaded from a file in the
following way:
```@repl oscar
C = positive_hull([1 0; 0 1])
save_cone(C, "C.cone")
CC = load_cone("C.cone")
collect(rays(CC))
```
The file is in JSON format and contains all previously gathered data belonging
to the underlying polymake object. In particular, this file can now be read by
both polymake and Oscar.

```@docs
save_cone(C::Cone, filename::String)
load_cone(filename::String)
```

## Auxiliary functions
```@docs
ambient_dim(C::Cone)
hilbert_basis(C::Cone)
codim(C::Cone)
dim(C::Cone)
ispointed(C::Cone)
isfulldimensional(C::Cone)
lineality_space(C::Cone)
nrays(C::Cone)
rays(C::Cone)
```
