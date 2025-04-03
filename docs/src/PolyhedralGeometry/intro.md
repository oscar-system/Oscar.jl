```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Introduction

The polyhedral geometry part of OSCAR provides functionality for handling
- convex polytopes, unbounded polyhedra and cones
- polyhedral fans
- linear and mixed integer programs

General textbooks offering details on theory and algorithms include:
- [JT13](@cite)
- [Sch86](@cite)
- [Zie95](@cite)


## Tutorials

We encourage you to take a look at the tutorials on polyhedral geometry in
OSCAR, which can be found [here](https://www.oscar-system.org/tutorials/PolyhedralGeometry/).


## Scalar types

The objects from polyhedral geometry operate on a given type, which (usually) resembles a field. This is indicated by the template parameter, e.g. the properties of a `Polyhedron{QQFieldElem}` are rational numbers of type `QQFieldElem`, if applicable.
Supported scalar types are `FieldElem` and `Float64`, but some functionality might not work properly if the parent `Field` does not satisfy certain mathematic conditions, like being ordered.
When constructing a polyhedral object from scratch, for the "simpler" types `QQFieldElem` and `Float64` it suffices to pass the `Type`, but more complex `FieldElem`s require a parent `Field` object. This can be set by either passing the desired `Field` instead of the type, or by inserting the type and have a matching `FieldElem` in your input data. If no type or field is given, the scalar type defaults to `QQFieldElem`.

The parent `Field` of the coefficients of an object `O` with coefficients of type `T` can be retrieved with the `coefficient_field` function, and it holds `elem_type(coefficient_field(O)) == T`.

```@docs
coefficient_field(x::PolyhedralObject)
```

!!! warning
    Support for fields other than the rational numbers is currently in an experimental stage.

These three lines result in the same polytope over rational numbers. Besides the general support mentioned above, naming a `Field` explicitly is encouraged because it allows user control and increases efficiency.
```jldoctest
julia> P = convex_hull(QQ, [1 0 0; 0 0 1]) # passing a `Field` always works
Polyhedron in ambient dimension 3

julia> P == convex_hull(QQFieldElem, [1 0 0; 0 0 1]) # passing the type works for `QQFieldElem` and `Float64` only
true

julia> P == convex_hull([1 0 0; 0 0 1]) # `Field` defaults to `QQ`
true

```

## Type compatibility

When working in polyhedral geometry it can prove advantageous to have various
input formats for the same kind of re-occurring quantitative input information.
This example shows three different ways to write the points whose convex hull
is to be computed, all resulting in identical `Polyhedron` objects:

```jldoctest
julia> P = convex_hull([1 0 0; 0 0 1])
Polyhedron in ambient dimension 3

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

There are two specialized `Vector`-like types, `PointVector` and `RayVector`, which commonly are returned by functions from Polyhedral Geometry. These can also be manually constructed:

```@docs
point_vector
ray_vector
```

While `RayVector`s can not be used do describe `PointVector`s (and vice versa),
matrices are generally allowed.

The primitive generator, also called minimal generator, of a ray can be
accessed as follows:

```@docs
primitive_generator(r::AbstractVector{T}) where T<:RationalUnion
primitive_generator_with_scaling_factor(r::AbstractVector{T}) where T<:RationalUnion
```

`AbstractCollection[PointVector]` can be given as:

Type                               | A `PointVector` corresponds to...
:--------------------------------- | :-------------------------------------------------------
`AbstractVector{<:PointVector}`    | an element of the vector.
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

Similar to points and rays, there are types `AffineHalfspace`, `LinearHalfspace`, `AffineHyperplane` and `LinearHyperplane`:

```@docs
affine_halfspace
linear_halfspace
affine_hyperplane
linear_hyperplane
```

These collections allow to mix up affine halfspaces/hyperplanes and their linear
counterparts, but note that an error will be produced when trying to convert
an affine description with bias not equal to zero to a linear description.

`AbstractCollection[LinearHalfspace]` can be given as:

Type                                   | A `LinearHalfspace` corresponds to...
:------------------------------------- | :----------------------------------------------------------
`AbstractVector{<:Halfspace}`          | an element of the vector.
`AbstractMatrix`/`MatElem` `A`         | the halfspace with normal vector `A[i, :]`.
`AbstractVector{<:AbstractVector}` `A` | the halfspace with normal vector `A[i]`.
`SubObjectIterator{<:Halfspace}`       | an element of the iterator.

`AbstractCollection[LinearHyperplane]` can be given as:

Type                                   | A `LinearHyperplane` corresponds to...
:------------------------------------- | :-----------------------------------------------------------
`AbstractVector{<:Hyperplane}`         | an element of the vector.
`AbstractMatrix`/`MatElem` `A`         | the hyperplane with normal vector `A[i, :]`.
`AbstractVector{<:AbstractVector}` `A` | the hyperplane with normal vector `A[i]`.
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


## `IncidenceMatrix`

Some methods will require input or return output in form of an `IncidenceMatrix`.

```@docs
IncidenceMatrix
```

The unique nature of the `IncidenceMatrix` allows for different ways of construction:

```@docs
incidence_matrix
```

From the examples it can be seen that this type supports `julia`'s matrix functionality. There are also functions to retrieve specific rows or columns as a `Set` over the non-zero indices.

```@docs
row(i::IncidenceMatrix, n::Int)
column(i::IncidenceMatrix, n::Int)
```

A typical application is the assignment of rays to the cones of a polyhedral
fan for its construction, see [`polyhedral_fan`](@ref).


## Visualization

Lower dimensional polyhedral objects can be visualized through polymake's backend.

```@docs
visualize(P::Union{Polyhedron{<:Union{Float64,FieldElem}}, Cone{<:Union{Float64,FieldElem}}, PolyhedralFan{<:Union{Float64,FieldElem}}, PolyhedralComplex{<:Union{Float64,FieldElem}}, SubdivisionOfPoints{<:Union{Float64,FieldElem}}, SimplicialComplex}; kwargs...)
visualize(::Vector)
```


## Serialization

Most objects from the polyhedral geometry section can be saved through the
polymake interface in the background. These functions are documented in the
subsections on the different objects. The format of the files is JSON and you
can find details of the specification
[here](https://polymake.org/schemas/data.json).

More details on the serialization, albeit concerning the older XML format, can be
found in [GHJ16](@cite). Even though the underlying format changed to JSON, the
abstract mathematical structure of the data files is still the same.


## Contact

Please direct questions about this part of OSCAR to the following people:
* [Taylor Brysiewicz](https://sites.google.com/view/taylorbrysiewicz/home),
* [Michael Joswig](https://page.math.tu-berlin.de/~joswig/),
* [Lars Kastner](https://lkastner.github.io/),
* Benjamin Lorenz.

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).
