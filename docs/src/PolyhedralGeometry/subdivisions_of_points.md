```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["subdivision_of_points.md"]
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
SubdivisionOfPoints(Points::Union{Oscar.MatElem,AbstractMatrix}, Incidence::IncidenceMatrix)
SubdivisionOfPoints(Points::Union{Oscar.MatElem,AbstractMatrix}, Weights::AbstractVector)
```

From a subdivision of points one can construct the
[`secondary_cone(SOP::SubdivisionOfPoints) `](@ref), i.e. the cone that is the
closure of the set of all weight vectors realizing that subdivision.

## Saving and loading

Objects of type `SubdivisionsOfPoints` can be saved to a file and loaded from a
file in the following way:
```jldoctest
julia> moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2]

julia> moaeincidence = IncidenceMatrix([[4,5,6],[1,4,2],[2,4,5],[2,3,5],[3,5,6],[1,3,6],[1,4,6]])

julia> MOAE = SubdivisionOfPoints(moaepts, moaeincidence)

julia> save("moae.sop", MOAE);

julia> SOP = load("moae.sop");

julia> is_regular(SOP)

```
The file is in JSON format and contains all previously gathered data belonging
to the underlying polymake object. In particular, this file can now be read by
both polymake and OSCAR.


## Auxiliary functions
```@docs
ambient_dim(SOP::SubdivisionOfPoints)
is_regular(SOP::SubdivisionOfPoints)
maximal_cells
min_weights
n_maximal_cells(SOP::SubdivisionOfPoints)
points(SOP::SubdivisionOfPoints)
secondary_cone
```
