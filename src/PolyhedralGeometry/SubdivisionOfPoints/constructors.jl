###############################################################################
###############################################################################
### Definition and constructors
###############################################################################
###############################################################################

struct SubdivisionOfPoints{T}
   pm_subdivision::Polymake.BigObject

   SubdivisionOfPoints{T}(pm::Polymake.BigObject) where T<:scalar_types = new{T}(pm)
end


# default scalar type: `QQFieldElem`
subdivision_of_points(points::AbstractCollection[PointVector], cells) =
  subdivision_of_points(QQFieldElem, points, cells)
subdivision_of_points(points::AbstractCollection[PointVector], weights::AbstractVector) =
  subdivision_of_points(QQFieldElem, points, weights)

# Automatic detection of corresponding OSCAR scalar type;
# Avoid, if possible, to increase type stability
subdivision_of_points(p::Polymake.BigObject) = SubdivisionOfPoints{detect_scalar_type(SubdivisionOfPoints, p)}(p)

@doc raw"""
    subdivision_of_points(points, cells)

# Arguments
- `points::AbstractCollection[PointVector]`: Points generating the cells of the
  subdivision; encoded row-wise as representative vectors.
- `cells::IncidenceMatrix`: An incidence matrix; there is a 1 at position (i,j)
  if cell i contains point j, and 0 otherwise.

A subdivision of points formed from points and cells made of these points. The
cells are given as an IncidenceMatrix, where the columns represent the points
and the rows represent the cells.

# Examples
The following is the famous "mother of all examples" (MOAE) non-regular
triangulation.
```jldoctest
julia> moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2];

julia> moaeimnonreg0 = IncidenceMatrix([[4,5,6],[1,4,2],[2,4,5],[2,3,5],[3,5,6],[1,3,6],[1,4,6]]);

julia> MOAE = subdivision_of_points(moaepts, moaeimnonreg0)
Subdivision of points in ambient dimension 3
```
"""
function subdivision_of_points(::Type{T}, points::AbstractCollection[PointVector], cells::IncidenceMatrix) where T<:scalar_types
   arr = @Polymake.convert_to Array{Set{Int}} Polymake.common.rows(cells)
   SubdivisionOfPoints{T}(Polymake.fan.SubdivisionOfPoints{_scalar_type_to_polymake(T)}(
      POINTS = homogenize(points,1),
      MAXIMAL_CELLS = arr,
   ))
end


@doc raw"""
    subdivision_of_points(points, weights)

# Arguments
- `points::AbstractCollection[PointVector]`: Points generating the cells of the
  subdivision; encoded row-wise as representative vectors.
- `weights::AbstractVector`: A vector with one entry for every point indicating
  the height of this point.

A subdivision of points formed by placing every point at the corresponding
height, then taking the convex hull and then only considering those cells
corresponding to faces visible from below ("lower envelope").

# Examples
We use the MOAE points, but give a weight vector instead of cells:
```jldoctest
julia> moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2];

julia> SOP = subdivision_of_points(moaepts, [1,1,1,1,1,1])
Subdivision of points in ambient dimension 3

julia> n_maximal_cells(SOP)
1
```
"""
function subdivision_of_points(::Type{T}, points::AbstractCollection[PointVector], weights::AbstractVector) where T<:scalar_types
   SubdivisionOfPoints{T}(Polymake.fan.SubdivisionOfPoints{_scalar_type_to_polymake(T)}(
      POINTS = homogenize(points,1),
      WEIGHTS = weights,
   ))
end

"""
    pm_object(SOP::SubdivisionOfPoints)

Get the underlying polymake object, which can be used via Polymake.jl.
"""
pm_object(SOP::SubdivisionOfPoints) = SOP.pm_subdivision


#Same construction for when the user provides maximal cells
function subdivision_of_points(::Type{T}, points::AbstractCollection[PointVector], cells::Vector{Vector{Int64}}) where T<:scalar_types
   subdivision_of_points(T, points, IncidenceMatrix(cells))
end

#Same construction for when the user gives Matrix{Bool} as incidence matrix
function subdivision_of_points(::Type{T}, Points::Union{Oscar.MatElem,AbstractMatrix}, cells::Matrix{Bool}) where T<:scalar_types
   subdivision_of_points(T, points, IncidenceMatrix(Polymake.IncidenceMatrix(cells)))
end

###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, SOP::SubdivisionOfPoints{T}) where T<:scalar_types
    print(io, "Subdivision of points in ambient dimension $(ambient_dim(SOP))")
    T != QQFieldElem && print(io, " with $T type coefficients")
end
