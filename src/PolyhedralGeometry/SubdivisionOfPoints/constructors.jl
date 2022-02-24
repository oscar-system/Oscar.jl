###############################################################################
###############################################################################
### Definition and constructors
###############################################################################
###############################################################################

struct SubdivisionOfPoints
   pm_subdivision::Polymake.BigObject
   function SubdivisionOfPoints(pm::Polymake.BigObject)
      return new(pm)
   end
end

@doc Markdown.doc"""
    SubdivisionOfPoints(points::Union{Oscar.MatElem,AbstractMatrix}, cells::IncidenceMatrix)

# Arguments
- `points::Matrix`: Points generating the cells of the subdivision; encoded row-wise as representative vectors.
- `cells::IncidenceMatrix`: An incidence matrix; there is a 1 at position (i,j) if cell i contains point j, and 0 otherwise.

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
function SubdivisionOfPoints(points::Union{Oscar.MatElem,AbstractMatrix}, cells::IncidenceMatrix)
   arr = @Polymake.convert_to Array{Set{Int}} Polymake.common.rows(cells)
   SubdivisionOfPoints(Polymake.fan.SubdivisionOfPoints{Polymake.Rational}(
      POINTS = homogenize(points,1),
      MAXIMAL_CELLS = arr,
   ))
end


@doc Markdown.doc"""
    SubdivisionOfPoints(points::Union{Oscar.MatElem,AbstractMatrix}, weights::AbstractVector)

# Arguments
- `points::Matrix`: Points generating the cells of the subdivision; encoded row-wise as representative vectors.
- `weights::AbstractVector`: A vector with one entry for every point indicating the height of this point.

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
function SubdivisionOfPoints(points::Union{Oscar.MatElem,AbstractMatrix}, weights::AbstractVector)
   SubdivisionOfPoints(Polymake.fan.SubdivisionOfPoints{Polymake.Rational}(
      POINTS = homogenize(points,1),
      WEIGHTS = weights,
   ))
end

"""
    pm_object(SOP::SubdivisionOfPoints)

Get the underlying polymake object, which can be used via Polymake.jl.
"""
pm_object(SOP::SubdivisionOfPoints) = SOP.pm_subdivision




#Same construction for when the user gives Matrix{Bool} as incidence matrix
function SubdivisionOfPoints(points::Union{Oscar.MatElem,AbstractMatrix}, cells::Matrix{Bool})
   SubdivisionOfPoints(points,IncidenceMatrix(Polymake.IncidenceMatrix(cells)))
end





#Same construction for when the user provides maximal cells
function SubdivisionOfPoints(points::Union{Oscar.MatElem,AbstractMatrix}, cells::Vector{Vector{Int64}})
   SubdivisionOfPoints(points,IncidenceMatrix(cells))
end

###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, SOP::SubdivisionOfPoints)
    print(io, "A subdivision of points in ambient dimension $(ambient_dim(SOP))")
end
