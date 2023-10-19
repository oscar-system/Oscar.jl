```@meta
CurrentModule = Oscar
```
# Double complexes -- the user's interface
We briefly review the mathematical notion of a double complex. 
Let ``\mathcal A`` be an Abelian category. A double complex 
``D_{\bullet, \bullet}`` consists of a collection of objects ``D_{i, j}`` in 
``\mathcal A`` with indices ``(i, j) \in \mathbb Z^2`` and usually arranged 
in a matrix-like grid, together with two collections 
of morphisms ``D_{i, j} \to D_{i, j \pm 1}``, the *horizontal* morphisms, and 
``D_{i, j} \to D_{i \pm 1, j}``, the *vertical* morphisms, so that both 
the rows and the columns of ``D_{\bullet, \bullet}`` are complexes in the 
classical sense and such that all resulting squares of maps commute.

In practice one usually encounters complexes which are *bounded* in the sense 
that outside some specified area of indices ``(i, j) \in \mathbb Z^2`` the entries 
``D_{i, j}`` are all zero. Such entries are then usually omitted. 

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
    horizontal_typ(dc::AbsDoubleComplexOfMorphisms)
    vertical_typ(dc::AbsDoubleComplexOfMorphisms)
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
In case any of these bounds does not exist, say `has_right_bound(D) == false` 
for some `AbsDoubleComplexOfMorphisms` `D`, then every request `D[i, j]` with ``i \gg 0`` arbitrarily 
big (but at least greater or equal to any potential `left_bound` of `D`) 
should be considered **to be legitimate**. Thus the programmer is responsible to indicate the 
bounds for the indices `(i, j)` if 
for their individual implementation of a double complex 
there are limitations for legitimate requests `D[i, j]`.

If any of the above bounds exist, asking for `D[i, j]` beyond those bounds might 
nevertheless be legitimate. Consider, 
for instance, the case where the rows of a double complex `D` are free resolutions of 
modules. Sometimes these are not computed at creation, but only on request of a specific entry 
`D[i, j]`. In this case the `right_bound(D)` indicates the supremum of indices `i` 
for which `D[i, j]` has already been computed, but it should nevertheless be allowed to also 
ask for the entry `D[i+1, j]`. Whether or not such requests are admissible can be checked using
```@docs
    extends_right(D::AbsDoubleComplexOfMorphisms)
    extends_left(D::AbsDoubleComplexOfMorphisms)
    extends_up(D::AbsDoubleComplexOfMorphisms)
    extends_down(D::AbsDoubleComplexOfMorphisms)
    is_complete(D::AbsDoubleComplexOfMorphisms)
    is_horizontally_complete(D::AbsDoubleComplexOfMorphisms)
    is_vertically_complete(D::AbsDoubleComplexOfMorphisms)
```

