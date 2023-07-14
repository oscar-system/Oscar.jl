```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

# Cones


## Introduction

Let $\mathbb{F}$ be an ordered field; the default is that
$\mathbb{F}=\mathbb{Q}$ is the field of rational numbers.

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
cone_from_inequalities
cone_from_equations
secondary_cone(SOP::SubdivisionOfPoints{T}) where T<:scalar_types
```

## Auxiliary functions
```@docs
ambient_dim(C::Cone)
Base.in(v::AbstractVector, C::Cone)
Base.issubset(C0::Cone{T}, C1::Cone{T}) where T<:scalar_types
f_vector(C::Cone)
hilbert_basis(C::Cone{QQFieldElem})
codim(C::Cone)
dim(C::Cone)
polarize(C::Cone{T}) where T<:scalar_types
intersect(C::Cone...)
is_pointed(C::Cone)
is_fulldimensional(C::Cone)
lineality_dim(C::Cone)
lineality_space(C::Cone{T}) where T<:scalar_types
nfacets(C::Cone)
nrays(C::Cone)
rays(C::Cone{T}) where T<:scalar_types
rays_modulo_lineality(C::Cone{T}) where T<:scalar_types
```
