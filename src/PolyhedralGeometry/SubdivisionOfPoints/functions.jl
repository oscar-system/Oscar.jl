###############################################################################
###############################################################################
### Functions taking a SubdivisionOfPoints
###############################################################################
###############################################################################

@doc Markdown.doc"""
    secondary_cone(SOP::SubdivisionOfPoints)

Return the secondary cone of a subdivision of points, the closure of all the
weight vectors inducing the given subdivision of points.

# Examples
For a non-regular subdivision, the secondary cone can still contain non-trivial
weights, but it will not be full-dimensional.
```jldoctest
julia> moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2];

julia> moaeimnonreg0 = IncidenceMatrix([[4,5,6],[1,4,2],[2,4,5],[2,3,5],[3,5,6],[1,3,6],[1,4,6]]);

julia> MOAE = SubdivisionOfPoints(moaepts, moaeimnonreg0)
A subdivision of points in ambient dimension 3

julia> C = secondary_cone(MOAE)
A polyhedral cone in ambient dimension 6

julia> dim(C)
4
```
"""
function secondary_cone(SOP::SubdivisionOfPoints{T}) where T<:scalar_types
   Cone{T}(Polymake.fan.secondary_cone(pm_object(SOP)))
end


