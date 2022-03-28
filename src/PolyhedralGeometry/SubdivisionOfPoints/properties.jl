###############################################################################
###############################################################################
### Iterators
###############################################################################
###############################################################################


@doc Markdown.doc"""
    points(SOP::SubdivisionOfPoints)

Return the points of the subdivision of points, `SOP`.

# Examples
Display the points of the "mother of all examples" non-regular triangulation.
```jldoctest
julia> moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2];

julia> moaeimnonreg0 = IncidenceMatrix([[4,5,6],[1,4,2],[2,4,5],[2,3,5],[3,5,6],[1,3,6],[1,4,6]]);

julia> MOAE = SubdivisionOfPoints(moaepts, moaeimnonreg0);

julia> points(MOAE)
6-element SubObjectIterator{PointVector{fmpq}}:
 [4, 0, 0]
 [0, 4, 0]
 [0, 0, 4]
 [2, 1, 1]
 [1, 2, 1]
 [1, 1, 2]
"""
function points(SOP::SubdivisionOfPoints)
    return SubObjectIterator{PointVector{fmpq}}(pm_object(SOP), _point, size(pm_object(SOP).POINTS, 1))
end

_point(::Type{PointVector{fmpq}}, SOP::Polymake.BigObject, i::Base.Integer) = PointVector{fmpq}(SOP.POINTS[i, 2:end])

_point_matrix(::Val{_point}, SOP::Polymake.BigObject; homogenized=false) = SOP.POINTS[:, (homogenized ? 1 : 2):end]

_matrix_for_polymake(::Val{_point}) = _point_matrix





@doc Markdown.doc"""
    maximal_cells(SOP::SubdivisionOfPoints)

Return an iterator over the maximal cells of `SOP`.

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

julia> moaeimnonreg0 = IncidenceMatrix([[4,5,6],[1,4,2],[2,4,5],[2,3,5],[3,5,6],[1,3,6],[1,4,6]])
7×6 IncidenceMatrix
[4, 5, 6]
[1, 2, 4]
[2, 4, 5]
[2, 3, 5]
[3, 5, 6]
[1, 3, 6]
[1, 4, 6]


julia> MOAE = SubdivisionOfPoints(moaepts, moaeimnonreg0);

julia> maximal_cells(MOAE)
7-element SubObjectIterator{Polymake.Set{Int64}}:
 pm::Set<long, pm::operations::cmp>
{4 5 6}
 pm::Set<long, pm::operations::cmp>
{1 2 4}
 pm::Set<long, pm::operations::cmp>
{2 4 5}
 pm::Set<long, pm::operations::cmp>
{2 3 5}
 pm::Set<long, pm::operations::cmp>
{3 5 6}
 pm::Set<long, pm::operations::cmp>
{1 3 6}
 pm::Set<long, pm::operations::cmp>
{1 4 6}
```
"""
function maximal_cells(SOP::SubdivisionOfPoints)
    return SubObjectIterator{Polymake.Set{Polymake.to_cxx_type(Int)}}(pm_object(SOP), _maximal_cell, size(pm_object(SOP).MAXIMAL_CELLS, 1))
end

_maximal_cell(::Type{Polymake.Set{Polymake.to_cxx_type(Int)}}, SOP::Polymake.BigObject, i::Base.Integer) = Polymake.row(SOP.MAXIMAL_CELLS, i)


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

julia> SOP = SubdivisionOfPoints(moaepts, [1,1,1,1,1,1]);

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

julia> SOP = SubdivisionOfPoints(moaepts, [1,1,1,1,1,1]);

julia> ambient_dim(SOP)
3
```
"""
ambient_dim(SOP::SubdivisionOfPoints) = pm_object(SOP).VECTOR_AMBIENT_DIM::Int - 1


@doc Markdown.doc"""
    npoints(SOP::SubdivisionOfPoints)

Return the number of points of a `SubdivisionOfPoints`.

# Examples
```jldoctest
julia> moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2];

julia> SOP = SubdivisionOfPoints(moaepts, [1,1,1,1,1,1]);

julia> npoints(SOP)
6
```
"""
npoints(SOP::SubdivisionOfPoints) = pm_object(SOP).N_POINTS::Int



###############################################################################
## Points properties
###############################################################################

@doc Markdown.doc"""
    min_weights(SOP::SubdivisionOfPoints)

Return the minimal weights inducing a subdivision of points. This method will
give an error if the input subdivision is non-regular.

# Examples
If all points have the same weight, then the 0-vector is minimal.
```jldoctest
julia> moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2];

julia> SOP = SubdivisionOfPoints(moaepts, [1,1,1,1,1,1]);

julia> min_weights(SOP)
pm::Vector<long>
0 0 0 0 0 0
```
"""
function min_weights(SOP::SubdivisionOfPoints)
   pm_object(SOP).MIN_WEIGHTS
end


@doc Markdown.doc"""
    maximal_cells_as_incidence_matrix(SOP::SubdivisionOfPoints)

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

julia> SOP = SubdivisionOfPoints(moaepts, [1,1,1,1,1,1])
A subdivision of points in ambient dimension 3

julia> maximal_cells_as_incidence_matrix(SOP)
1×6 IncidenceMatrix
[1, 2, 3, 4, 5, 6]
```
"""
function maximal_cells_as_incidence_matrix(SOP::SubdivisionOfPoints)
   IncidenceMatrix(pm_object(SOP).MAXIMAL_CELLS)
end

###############################################################################
## Boolean properties
###############################################################################

"""
    isregular(SOP::SubdivisionOfPoints)

Determine whether `SOP` is regular, i.e. can be given via a height function.

# Examples
This is the so-called "mother of all examples", a very famous non-regular
triangulation of six points.
```jldoctest
julia> moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2];

julia> moaeimnonreg0 = IncidenceMatrix([[4,5,6],[1,4,2],[2,4,5],[2,3,5],[3,5,6],[1,3,6],[1,4,6]]);

julia> MOAE = SubdivisionOfPoints(moaepts, moaeimnonreg0);

julia> isregular(MOAE)
false

julia> SOP = SubdivisionOfPoints(moaepts, [1,1,1,1,1,1]);

julia> isregular(SOP)
true
```
"""
isregular(SOP::SubdivisionOfPoints) = pm_object(SOP).REGULAR::Bool
