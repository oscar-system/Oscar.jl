```@meta
CurrentModule = Oscar
```
# Double complexes -- the user's interface
We briefly review the mathematical notion of a double complex. 
Let ``\mathcal A`` be an Abelian category. A double complex 
``D_{\bullet, \bullet}`` consists of a collection of objects ``D_{i, j}`` in 
``\mathcal A`` with indices ``(i, j) \in \mathbb Z^2`` and usually arranged 
in a matrix-like grid, together with two collections 
of morphisms ``D_{i, j} \to D_{i \pm 1, j}``, the *horizontal* morphisms, and 
``D_{i, j} \to D_{i, j \pm 1}``, the *vertical* morphisms, so that both 
the rows and the columns of ``D_{\bullet, \bullet}`` are complexes in the 
classical sense and such that all resulting squares of maps commute.

In practice one usually encounters complexes which are *bounded* in the sense 
that outside some specified area of indices ``(i, j) \in \mathbb Z^2`` the entries 
``D_{i, j}`` are all zero. Such entries are then usually omitted. 

## Basic getters and attributes
In OSCAR the generic functionality for double complexes is declared for the 
abstract type `AbsDoubleComplexOfMorphisms`. These functions comprise
```julia
  getindex(D::AbsDoubleComplexOfMorphisms, i::Int, j::Int) # Get the `(i,j)`-th entry of `D`
```
```@docs
    horizontal_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
    vertical_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
```
In which direction the maps in the rows and columns go can be asked with the following methods:
```@docs
    horizontal_direction(dc::AbsDoubleComplexOfMorphisms)
    vertical_direction(dc::AbsDoubleComplexOfMorphisms)
```
Double complexes can be bounded or unbounded. But even in the bounded case where only finitely 
many ``D_{i, j}`` can be non-zero, it is not clear that a concrete bound for the 
range of indices ``(i, j)`` of possibly non-zero entries is known. Then again, even in that case 
such a bound might be rather theoretical and it is still encouraged that double complexes be implemented 
lazy and not fill out the full grid of morphisms on construction. Thus a merely theoretical bound 
on the range of indices does not seem to be practically relevant. 

To accomodate all of these cases, we are forced to employ a rather general pattern for the 
admissible ranges of a double complex.
```@docs
    has_upper_bound(D::AbsDoubleComplexOfMorphisms)
    has_lower_bound(D::AbsDoubleComplexOfMorphisms)
    has_right_bound(D::AbsDoubleComplexOfMorphisms)
    has_left_bound(D::AbsDoubleComplexOfMorphisms)
```
If they exist, these bounds can be asked for using 
```@docs
    right_bound(D::AbsDoubleComplexOfMorphisms)
    left_bound(D::AbsDoubleComplexOfMorphisms)
    upper_bound(D::AbsDoubleComplexOfMorphisms)
    lower_bound(D::AbsDoubleComplexOfMorphisms)
```
The existence of the above bounds does **not** indicate whether the entry 
`D[i, j]` or the maps leaving from it have been or can be computed! 
In principal it is advised to check for this using 
```@docs 
    has_index(D::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
    can_compute_index(D::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
    has_horizontal_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
    can_compute_horizontal_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
    has_vertical_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
    can_compute_vertical_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
```
It is also possible to query whether or not a double complex 
is already complete in the sense that it knows about all of its 
non-zero entries.
```@docs
    is_complete(D::AbsDoubleComplexOfMorphisms)
```

## Generic functionality
```@docs
    total_complex(D::AbsDoubleComplexOfMorphisms)
```

## Constructors
```@docs 
    tensor_product(C1::ComplexOfMorphisms{ChainType}, C2::ComplexOfMorphisms{ChainType}) where {ChainType}
```

