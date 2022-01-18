###############################################################################
###############################################################################
### Definition and constructors
###############################################################################
###############################################################################

struct PolyhedralComplex
    pm_complex::Polymake.BigObject
    function PolyhedralComplex(pm::Polymake.BigObject)
        return new(pm)
    end
end

pm_object(pc::PolyhedralComplex) = pc.pm_complex


@doc Markdown.doc"""
    PolyhedralComplex(polyhedra, vr, far_vertices, L)

# Arguments
- `polyhedra::IncidenceMatrix`: An incidence matrix; there is a 1 at position
  (i,j) if the ith polytope contains point j and 0 otherwise.
- `vr::Matrix`: The points whose convex hulls make up the polyhedral
  complex. This matrix also contains the far vertices.
- `far_vertices::Vector{Int}`: Vector containing the indices of the rows
  corresponding to the far vertices in `vr`.
- `L::Matrix`: Generators of the lineality space of the polyhedral complex.

A polyhedral complex formed from points, rays, and lineality combined into
polyhedra indicated by an incidence matrix, where the columns represent the
points and the rows represent the polyhedra.

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2,3],[1,3,4]])
2×4 IncidenceMatrix
[1, 2, 3]
[1, 3, 4]


julia> vr = [0 0; 1 0; 1 1; 0 1]
4×2 Matrix{Int64}:
 0  0
 1  0
 1  1
 0  1

julia> PC = PolyhedralComplex(IM, vr)
A polyhedral complex in ambient dimension 2
```
"""
function PolyhedralComplex(
                polyhedra::IncidenceMatrix, 
                vr::Union{SubObjectIterator{<:Union{PointVector,PointVector}}, Oscar.MatElem, AbstractMatrix}, 
                far_vertices::Union{Vector{Int}, Nothing} = nothing, 
                L::Union{SubObjectIterator{<:RayVector}, 
                Oscar.MatElem, AbstractMatrix, Nothing} = nothing
            )
    LM = isnothing(L) || isempty(L) ? Polymake.Matrix{Polymake.Rational}(undef, 0, size(vr, 2)) : L

    # Rays and Points are homogenized and combined and
    points = homogenize(vr, 1)
    # If some vertices are far vertices, give them a leading 0
    if !isnothing(far_vertices)
        points[far_vertices,1] .= 0
    end

    # Lineality is homogenized
    lineality = homogenize(LM, 0)

    PolyhedralComplex(Polymake.fan.PolyhedralComplex{Polymake.Rational}(
        POINTS = points,
        INPUT_LINEALITY = lineality,
        INPUT_CONES = polyhedra,
    ))
end

# This makes sure that PolyhedralComplex(maximal_polyhedra(PC)) returns an Oscar PolyhedralComplex,
PolyhedralComplex(iter::SubObjectIterator{Polyhedron}) = PolyhedralComplex(iter.Obj)

###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, PC::PolyhedralComplex)
    try
        ad = ambient_dim(PC)
        print(io, "A polyhedral complex in ambient dimension $(ad)")
    catch e
        print(io, "A polyhedral complex without ambient dimension")
    end
end
