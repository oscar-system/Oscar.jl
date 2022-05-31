```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["intro.md"]
```

# Introduction

The polyhedral geometry part of OSCAR provides functionality for handling
- convex polytopes, unbounded polyhedra and cones
- polyhedral fans
- linear programs

General textbooks offering details on theory and algorithms include:
- [JT13](@cite)
- [Zie95](@cite)


## Type compatibility

When working in polyhedral geometry it can prove advantageous to have various
input formats for the same kind of re-occuring quantitative input information.
This example shows three different ways to name the points whose convex hull
is to be computed, all resulting in identical `Polyhedron` objects:

```
julia> P = convex_hull([1 0 0; 0 0 1])
A polyhedron in ambient dimension 3

julia> P == convex_hull([[1, 0, 0], [0, 0, 1]])
true

julia> P == convex_hull(vertices(P))
true
```

`convex_hull` is only one of many functions and constructors supporting this
behavior, and there are also more types that can be described this way besides
`PointVector`. Whenever the docs state an argument is required to be of type
`AbstractCollection[ElType]` (where `ElType` is the `Oscar` type of single
instances described in this collection), the user can choose the input to follow
any of the corresponding notions below.


### Vectors

While `RayVector`s can not be used do describe `PointVector`s (and vice versa),
matrices are generally allowed.

`AbstractCollection[PointVector]` can be given as:

Type                               | A `PointVector` corresponds to...
:--------------------------------- | :-------------------------------------------------------
`AbstractVector{<:ṔointVector}`    | an element of the vector.
`AbstractVector{<:AbstractVector}` | an element of the vector.
`AbstractMatrix`/`MatElem`         | a row of the matrix.
`AbstractVector`/`PointVector`     | the vector itself (only one `PointVector` is described).
`SubObjectIterator{<:PointVector}` | an element of the iterator.

`AbstractCollection[RayVector]` can be given as:

Type                               | A `RayVector` corresponds to...
:--------------------------------- | :-----------------------------------------------------
`AbstractVector{<:RayVector}`      | an element of the vector.
`AbstractVector{<:AbstractVector}` | an element of the vector.
`AbstractMatrix`/`MatElem`         | a row of the matrix.
`AbstractVector`/`RayVector`       | the vector itself (only one `RayVector` is described).
`SubObjectIterator{<:RayVector}`   | an element of the iterator.


### Halfspaces and Hyperplanes

These collections allow to mix up affine halfspaces/hyperplanes and their linear
counterparts, but note that an error will be produced when trying to convert
an affine description with bias not equal to zero to a linear description.

`AbstractCollection[LinearHalfspace]` can be given as:

Type                                   | A `LinearHalfspace` corresponds to...
:------------------------------------- | :----------------------------------------------------------
`AbstractVector{<:Halfspace}`          | an element of the vector.
`AbstractMatrix`/`MatElem` `A`         | the halfspace with normal vector `A[i, :]`.
`SubObjectIterator{<:Halfspace}`       | an element of the iterator.

`AbstractCollection[LinearHyperplane]` can be given as:

Type                                   | A `LinearHyperplane` corresponds to...
:------------------------------------- | :-----------------------------------------------------------
`AbstractVector{<:Hyperplane}`         | an element of the vector.
`AbstractMatrix`/`MatElem` `A`         | the hyperplane with normal vector `A[i, :]`.
`SubObjectIterator{<:Hyperplane}`      | an element of the iterator.

`AbstractCollection[AffineHalfspace]` can be given as:

Type                                   | An `AffineHalfspace` corresponds to...
:------------------------------------- | :-----------------------------------------------------------------
`AbstractVector{<:Halfspace}`          | an element of the vector.
`Tuple` over matrix `A` and vector `b` | the affine halfspace with normal vector `A[i, :]` and bias `b[i]`.
`SubObjectIterator{<:Halfspace}`       | an element of the iterator.

`AbstractCollection[AffineHyperplane]` can be given as:

Type                                   | An `AffineHyperplane` corresponds to...
:------------------------------------- | :------------------------------------------------------------------
`AbstractVector{<:Hyperplane}`         | an element of the vector.
`Tuple` over matrix `A` and vector `b` | the affine hyperplane with normal vector `A[i, :]` and bias `b[i]`.
`SubObjectIterator{<:Hyperplane}`      | an element of the iterator.


## Serialization

Most objects from the polyhedral geometry section can be saved through the
polymake interface in the background. These functions are documented in the
subsections on the different objects. The format of the files is JSON and you
can find details of the specification
[here](https://polymake.org/schemas/data.json).

More details on the serialization, albeit concerning the older XML format, can be
found in [GHJ16](@cite). Even though the underlying format changed to JSON, the
abstract mathematical structure of the data files is still the same.
