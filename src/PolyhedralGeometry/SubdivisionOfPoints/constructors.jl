###############################################################################
###############################################################################
### Definition and constructors
###############################################################################
###############################################################################

struct SubdivisionOfPoints{T}
   pm_subdivision::Polymake.BigObject
   
   SubdivisionOfPoints{T}(pm::Polymake.BigObject) where T<:scalar_types = new{T}(pm)
end

# default scalar type: `fmpq`
SubdivisionOfPoints(x...) = SubdivisionOfPoints{fmpq}(x...)

# Automatic detection of corresponding OSCAR scalar type;
# Avoid, if possible, to increase type stability
SubdivisionOfPoints(p::Polymake.BigObject) = SubdivisionOfPoints{detect_scalar_type(SubdivisionOfPoints, p)}(p)

@doc Markdown.doc"""
    SubdivisionOfPoints(Points, Cells)

# Arguments
- `Points::Matrix`: Points generating the cells of the subdivision; encoded row-wise as representative vectors.
- `Cells::IncidenceMatrix`: An incidence matrix; there is a 1 at position (i,j) if cell i contains point j, and 0 otherwise.

A subdivision of points formed from points and cells made of these points. The
cells are given as an IncidenceMatrix, where the columns represent the points
and the rows represent the cells.

# Examples
The following is the famous "mother of all examples" (moae) non-regular
triangulation.
```jldoctest
julia> moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2];

julia> moaeimnonreg0 = IncidenceMatrix([[4,5,6],[1,4,2],[2,4,5],[2,3,5],[3,5,6],[1,3,6],[1,4,6]]);

julia> MOAE = SubdivisionOfPoints(moaepts, moaeimnonreg0)
A subdivision of points in ambient dimension 3
```
"""
function SubdivisionOfPoints{T}(Points::Union{Oscar.MatElem,AbstractMatrix}, Incidence::IncidenceMatrix) where T<:scalar_types
   arr = @Polymake.convert_to Array{Set{Int}} Polymake.common.rows(Incidence)
   SubdivisionOfPoints{T}(Polymake.fan.SubdivisionOfPoints{scalar_type_to_polymake[T]}(
      POINTS = homogenize(Points,1),
      MAXIMAL_CELLS = arr,
   ))
end


@doc Markdown.doc"""
    SubdivisionOfPoints(Points, Weights)

# Arguments
- `Points::Matrix`: Points generating the cells of the subdivision; encoded row-wise as representative vectors.
- `Weights::AbstractVector`: A vector with one entry for every point indicating the height of this point.

A subdivision of points formed by placing every point at the corresponding
height, then taking the convex hull and then only considering those cells
corresponding to faces visible from below ("lower envelope").

# Examples
We use the MOAE points, but give a weight vector instead of cells:
```jldoctest
julia> moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2];

julia> SOP = SubdivisionOfPoints(moaepts, [1,1,1,1,1,1])
A subdivision of points in ambient dimension 3

julia> n_maximal_cells(SOP)
1
```
"""
function SubdivisionOfPoints{T}(Points::Union{Oscar.MatElem,AbstractMatrix}, Weights::AbstractVector) where T<:scalar_types
   SubdivisionOfPoints{T}(Polymake.fan.SubdivisionOfPoints{scalar_type_to_polymake[T]}(
      POINTS = homogenize(Points,1),
      WEIGHTS = Weights,
   ))
end

"""
    pm_object(SOP::SubdivisionOfPoints)

Get the underlying polymake object, which can be used via Polymake.jl.
"""
pm_object(SOP::SubdivisionOfPoints) = SOP.pm_subdivision




#Same construction for when the user gives Matrix{Bool} as incidence matrix
function SubdivisionOfPoints{T}(Points::Union{Oscar.MatElem,AbstractMatrix}, Incidence::Matrix{Bool}) where T<:scalar_types
   SubdivisionOfPoints{T}(Points,IncidenceMatrix(Polymake.IncidenceMatrix(Incidence)))
end

###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, SOP::SubdivisionOfPoints{T}) where T<:scalar_types
    print(io, "A subdivision of points in ambient dimension $(ambient_dim(SOP))")
    T != fmpq && print(io, " with $T type coefficients")
end
