```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
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
Double complexes can be bounded or unbounded. It is important to note that even if such 
bounds exist and are known, this is a priori **not** related to whether or not 
certain entries are computable! I.e. even in the case of a bounded complex `dc` 
it might still be valid to call `dc[i, j]` beyond that bound. In general, one should 
use the following functions to determine whether or not it is legitimate to ask for a 
specific entry.
```@docs 
    has_index(D::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
    can_compute_index(D::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
    has_horizontal_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
    can_compute_horizontal_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
    has_vertical_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
    can_compute_vertical_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
```
Explicitly known bounds for the non-zero entries of a complex are nevertheless relevant for 
various generic functionalities. 
For example, computing a total complex is only possible in practice if one has an a priori estimate 
where the non-zero entries are located. For such purposes, we provide the following functionality:
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

