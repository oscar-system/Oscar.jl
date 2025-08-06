###############################################################################
###############################################################################
### Iterators
###############################################################################
###############################################################################

@doc raw"""
    points(SOP::SubdivisionOfPoints)

Return the points of the subdivision of points, `SOP`.

# Examples
Display the points of the "mother of all examples" non-regular triangulation.
```jldoctest
julia> moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2];

julia> moaeimnonreg0 = incidence_matrix([[4,5,6],[1,4,2],[2,4,5],[2,3,5],[3,5,6],[1,3,6],[1,4,6]]);

julia> MOAE = subdivision_of_points(moaepts, moaeimnonreg0);

julia> points(MOAE)
6-element SubObjectIterator{PointVector{QQFieldElem}}:
 [4, 0, 0]
 [0, 4, 0]
 [0, 0, 4]
 [2, 1, 1]
 [1, 2, 1]
 [1, 1, 2]
```
"""
function points(SOP::SubdivisionOfPoints{T}) where {T<:scalar_types}
  return SubObjectIterator{PointVector{T}}(SOP, _point, size(pm_object(SOP).POINTS, 1))
end

_point(
  U::Type{PointVector{T}}, SOP::SubdivisionOfPoints{T}, i::Base.Integer
) where {T<:scalar_types} =
  point_vector(coefficient_field(SOP), pm_object(SOP).POINTS[i, 2:end])::U

_point_matrix(::Val{_point}, SOP::SubdivisionOfPoints; homogenized=false) =
  pm_object(SOP).POINTS[:, (homogenized ? 1 : 2):end]

_matrix_for_polymake(::Val{_point}) = _point_matrix

@doc raw"""
    maximal_cells(SOP::SubdivisionOfPoints)

Return an iterator over the maximal cells of `SOP`.

Optionally `IncidenceMatrix` can be passed as a first argument to return the
incidence matrix specifying the maximal cells of `SOP`.

# Examples
Display the cells of the "mother of all examples" non-regular triangulation.
```jldoctest
julia> moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2]
6×3 Matrix{Int64}:
 4  0  0
 0  4  0
 0  0  4
 2  1  1
 1  2  1
 1  1  2

julia> moaeimnonreg0 = incidence_matrix([[4,5,6],[1,4,2],[2,4,5],[2,3,5],[3,5,6],[1,3,6],[1,4,6]])
7×6 IncidenceMatrix
 [4, 5, 6]
 [1, 2, 4]
 [2, 4, 5]
 [2, 3, 5]
 [3, 5, 6]
 [1, 3, 6]
 [1, 4, 6]

julia> MOAE = subdivision_of_points(moaepts, moaeimnonreg0);

julia> maximal_cells(MOAE)
7-element SubObjectIterator{Vector{Int64}}:
 [4, 5, 6]
 [1, 2, 4]
 [2, 4, 5]
 [2, 3, 5]
 [3, 5, 6]
 [1, 3, 6]
 [1, 4, 6]

julia> maximal_cells(IncidenceMatrix, MOAE)
7×6 IncidenceMatrix
 [4, 5, 6]
 [1, 2, 4]
 [2, 4, 5]
 [2, 3, 5]
 [3, 5, 6]
 [1, 3, 6]
 [1, 4, 6]
```
"""
maximal_cells(SOP::SubdivisionOfPoints) = maximal_cells(Vector{Int}, SOP)
function maximal_cells(::Type{Vector{Int}}, SOP::SubdivisionOfPoints)
  return SubObjectIterator{Vector{Int}}(
    SOP, _maximal_cell, size(pm_object(SOP).MAXIMAL_CELLS, 1)
  )
end

_maximal_cell(::Type{Vector{Int}}, SOP::SubdivisionOfPoints, i::Base.Integer) =
  Vector{Int}(Polymake.row(pm_object(SOP).MAXIMAL_CELLS, i))

###############################################################################
###############################################################################
### Access properties
###############################################################################
###############################################################################

###############################################################################
## Scalar properties
###############################################################################

"""
    n_maximal_cells(SOP::SubdivisionOfPoints)

Return the number of maximal cells of `SOP`.
# Examples
If all points have the same weight, there is only one cell.
```jldoctest
julia> moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2];

julia> SOP = subdivision_of_points(moaepts, [1,1,1,1,1,1]);

julia> n_maximal_cells(SOP)
1
```
"""
n_maximal_cells(SOP::SubdivisionOfPoints) = pm_object(SOP).N_MAXIMAL_CELLS

"""
    ambient_dim(SOP::SubdivisionOfPoints)

Return the ambient dimension `SOP`, which is the dimension of the embedding
space.

# Examples
The ambient dimension of the MOAE is 3, independent of the subdivision chosen.
```jldoctest
julia> moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2];

julia> SOP = subdivision_of_points(moaepts, [1,1,1,1,1,1]);

julia> ambient_dim(SOP)
3
```
"""
ambient_dim(SOP::SubdivisionOfPoints) = pm_object(SOP).VECTOR_AMBIENT_DIM::Int - 1

@doc raw"""
    n_points(SOP::SubdivisionOfPoints)

Return the number of points of a `SubdivisionOfPoints`.

# Examples
```jldoctest
julia> moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2];

julia> SOP = subdivision_of_points(moaepts, [1,1,1,1,1,1]);

julia> n_points(SOP)
6
```
"""
n_points(SOP::SubdivisionOfPoints) = pm_object(SOP).N_POINTS::Int

###############################################################################
## Points properties
###############################################################################

@doc raw"""
    min_weights(SOP::SubdivisionOfPoints)

Return the minimal weights inducing a subdivision of points. This method will
give an error if the input subdivision is non-regular.

# Examples
If all points have the same weight, then the 0-vector is minimal.
```jldoctest
julia> moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2];

julia> SOP = subdivision_of_points(moaepts, [1,1,1,1,1,1]);

julia> min_weights(SOP)
6-element Vector{Int64}:
 0
 0
 0
 0
 0
 0
```
"""
min_weights(SOP::SubdivisionOfPoints) = Vector{Int}(pm_object(SOP).MIN_WEIGHTS)

@doc raw"""
    maximal_cells(IncidenceMatrix, SOP::SubdivisionOfPoints)

Return the maximal cells of `SOP` as an incidence matrix.

The rows of the output correspond to the maximal cells and the columns
correspond to the cells.

# Examples
If we give all points the same weight there is only one cell containing all
points.
```jldoctest
julia> moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2]
6×3 Matrix{Int64}:
 4  0  0
 0  4  0
 0  0  4
 2  1  1
 1  2  1
 1  1  2

julia> SOP = subdivision_of_points(moaepts, [1,1,1,1,1,1])
Subdivision of points in ambient dimension 3

julia> maximal_cells(IncidenceMatrix, SOP)
1×6 IncidenceMatrix
 [1, 2, 3, 4, 5, 6]
```
"""
maximal_cells(::Type{IncidenceMatrix}, SOP::SubdivisionOfPoints) =
  pm_object(SOP).MAXIMAL_CELLS

###############################################################################
## Boolean properties
###############################################################################

"""
    is_regular(SOP::SubdivisionOfPoints)

Determine whether `SOP` is regular, i.e. can be given via a height function.

# Examples
This is the so-called "mother of all examples", a very famous non-regular
triangulation of six points.
```jldoctest
julia> moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2];

julia> moaeimnonreg0 = incidence_matrix([[4,5,6],[1,4,2],[2,4,5],[2,3,5],[3,5,6],[1,3,6],[1,4,6]]);

julia> MOAE = subdivision_of_points(moaepts, moaeimnonreg0);

julia> is_regular(MOAE)
false

julia> SOP = subdivision_of_points(moaepts, [1,1,1,1,1,1]);

julia> is_regular(SOP)
true
```
"""
is_regular(SOP::SubdivisionOfPoints) = pm_object(SOP).REGULAR::Bool
