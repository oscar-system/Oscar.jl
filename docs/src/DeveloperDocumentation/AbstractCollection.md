# `AbstractCollection`

Allowing the user to pass input using several formats usually is handled within
`julia` by defining specialized methods for each function and argument type(s).
This can prove to be inefficient when the amount possible combinations of these
increases. `AbstractCollection` is a `Dict` meant to enable the user to profit
from a fixed interpretation describing different collections of mathematical
objects, while also simplifying the life of the developer, also resulting
in less code duplication.

## Idea
Commonly the same kind of information, e.g. an amount of `PointVector`s, is
accepted as argument for many different functions. The user can chose from
different types (coming with an interpretation of their content) to use when
calling one of these functions to describe the data. This data is then converted
to a type and format `Polymake.jl` (and thus indirectly the `polymake` kernel)
supports.

## Example
Usually in `polymake`, a collection of points is displayed as a matrix
of row-vectors. Such a matrix is always created from the input information. When
writing a new function accepting an object `x` of type
`AbstractCollection[PointVector]` (note that, with `AbstractCollection` being a
`Dict`, its entries are accessed using square brackets; the keys are the `Oscar`
types of the elements of the collection), the necessary conversion can (and
should) be called at the beginning. These conversion functions already exist and
support all of the types stated in [Type compatibility](@ref). In this case the
function is `homogenized_matrix(x, 1)`.

`RayVector`s and their collections work about the same; the main difference for
the programmer is that `homogenized_matrix(x, 0)` is called.

When looking at the beginning of the `convex_hull` method, the corresponding
conversions of the three arguments `V`, `R` and `L` can be seen:

```julia
function convex_hull(::Type{T}, V::AbstractCollection[PointVector], R::Union{AbstractCollection[RayVector], Nothing} = nothing, L::Union{AbstractCollection[RayVector], Nothing} = nothing; non_redundant::Bool = false) where T<:scalar_types
    # Rays and Points are homogenized and combined and
    # Lineality is homogenized
    points = stack(homogenized_matrix(V, 1), homogenized_matrix(R, 0))
    lineality = isnothing(L) || isempty(L) ? zero_matrix(QQ, 0, size(points,2)) : homogenized_matrix(L, 0)

    ...
end
```

## Conversion functions

So effectively supporting `AbstractCollection`s only requires to know when to
apply which conversion function. The following table explains this for
`AbstractCollection[T]`:

`T`                                  | Target format                               | Conversion function
:----------------------------------- | :------------------------------------------ | :------------------------------
`PointVector`                        | matrix of row-vectors                       | `homogenized_matrix(*, 1)`
`RayVector`                          | matrix of row-vectors (linear setting)      | `unhomogenized_matrix(*)`
`RayVector`                          | matrix of row-vectors (affine setting)      | `homogenized_matrix(*, 0)`
`LinearHalfspace`/`LinearHyperplane` | inequality/equation matrix (linear setting) | `linear_matrix_for_polymake(*)`
`AffineHalfspace`/`AffineHyperplane` | inequality/equation matrix (affine setting) | `affine_matrix_for_polymake(*)`
