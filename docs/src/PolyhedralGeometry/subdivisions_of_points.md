```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

# Subdivisions of Points

## Introduction

Let $\mathbb{F}$ be an ordered field; the default is that
$\mathbb{F}=\mathbb{Q}$ is the field of rational numbers and other fields are
not yet supported everywhere in the implementation.

A *subdivision of points* consists of

- a finite set $\mathcal{P}\subseteq\mathbb{F}^n$ of points; and
- a finite set of cells $\mathcal{S}\subseteq 2^{\mathcal{P}}$.

The cells are only allowed to intersect in common faces. In contrast to the
maximal cones of a polyhedral fan or the maximal polytopes of a polyhedral
complex, cells are allowed to have interior points here, i.e. they are not
given in terms of their vertices.


## Construction

There are two ways to construct a subdivision of points. First, one can specify
the cells directly. Second, one can assign a weight or height to every point,
take the convex hull and take the cells corresponding to facets visible from
below ("lower envelope"). Not every subdivision of points comes from a weight
vector, but if it does, it is called `regular`.


```@docs
subdivision_of_points(::Oscar.scalar_type_or_field, Points::Union{Oscar.MatElem,AbstractMatrix}, Incidence::IncidenceMatrix)
subdivision_of_points(::Oscar.scalar_type_or_field, Points::Union{Oscar.MatElem,AbstractMatrix}, Weights::AbstractVector)
```

From a subdivision of points one can construct the
[`secondary_cone(SOP::SubdivisionOfPoints) `](@ref), i.e. the cone that is the
closure of the set of all weight vectors realizing that subdivision.


## Auxiliary functions
```@docs
ambient_dim(SOP::SubdivisionOfPoints)
is_regular(SOP::SubdivisionOfPoints)
maximal_cells
min_weights
number_of_maximal_cells(SOP::SubdivisionOfPoints)
points(SOP::SubdivisionOfPoints{T}) where T<:scalar_types
secondary_cone
```
