###############################################################################
###############################################################################
### Iterators
###############################################################################
###############################################################################

#Should these be coded as SubObjectIterators?
#  For example, so we can apply point_matrix(points(SOP))?
@doc Markdown.doc"""
    points(SOP::SubdivisionOfPoints)

Return an iterator over the points of `SOP`.
"""
function points(SOP::SubdivisionOfPoints)
   PointSOPIterator(SOP)
end

struct PointSOPIterator
    SOP::SubdivisionOfPoints
end

function Base.iterate(iter::PointSOPIterator, index = 1)
    n_points = npoints(iter.SOP)
    if index > n_points
        return nothing
    end
    #dehomogenize doesn't check if the first entry is 0. I think this okay here.
    current_point = dehomogenize(pm_object(iter.SOP).POINTS[index,:])
    return (current_point, index + 1)
end
Base.length(iter::PointSOPIterator) = npoints(iter.SOP)




"""
    maximal_cells(SOP::SubdivisionOfPoints)

Return the maximal cells of `SOP`.

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
7×6 Matrix{Bool}:
 0  0  0  1  1  1
 1  1  0  1  0  0
 0  1  0  1  1  0
 0  1  1  0  1  0
 0  0  1  0  1  1
 1  0  1  0  0  1
 1  0  0  1  0  1

julia> MOAE = SubdivisionOfPoints(moaepts, moaeimnonreg0);

julia> for c in maximal_cells(MOAE)
       println(c)
       end
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


@doc Markdown.doc"""
    maximal_cells(SOP::SubdivisionOfPoints)

Return an iterator over the maximal cells of `SOP`.
"""
function maximal_cells(SOP::SubdivisionOfPoints)
   MaximalCellIterator(SOP)
end

struct MaximalCellIterator
    SOP::SubdivisionOfPoints
end

function Base.iterate(iter::MaximalCellIterator, index = 1)
    n_max_cells = n_maximal_cells(iter.SOP)
    if index > n_max_cells
        return nothing
    end
    current_cell = Polymake.row(pm_object(iter.SOP).MAXIMAL_CELLS, index)
    return (current_cell, index + 1)
end
Base.length(iter::MaximalCellIterator) = n_maximal_cells(iter.SOP)

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
