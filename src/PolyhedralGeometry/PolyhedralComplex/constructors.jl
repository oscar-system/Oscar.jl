###############################################################################
###############################################################################
### Definition and constructors
###############################################################################
###############################################################################

struct PolyhedralComplex{T}
    pm_complex::Polymake.BigObject
    
     PolyhedralComplex{T}(pm::Polymake.BigObject) where T<:scalar_types = new{T}(pm)
end

# default scalar type: `fmpq`
PolyhedralComplex(x...) = PolyhedralComplex{fmpq}(x...)

PolyhedralComplex(p::Polymake.BigObject) = PolyhedralComplex{detect_scalar_type(PolyhedralComplex, p)}(p)

pm_object(pc::PolyhedralComplex) = pc.pm_complex


@doc Markdown.doc"""
    PolyhedralComplex{T}(polyhedra, vr, far_vertices, L) where T<:scalar_types

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
function PolyhedralComplex{T}(
                polyhedra::IncidenceMatrix, 
                vr::Union{SubObjectIterator{<:Union{PointVector,PointVector}}, Oscar.MatElem, AbstractMatrix}, 
                far_vertices::Union{Vector{Int}, Nothing} = nothing, 
                L::Union{SubObjectIterator{<:RayVector}, 
                Oscar.MatElem, AbstractMatrix, Nothing} = nothing
            ) where T<:scalar_types
    LM = isnothing(L) || isempty(L) ? Polymake.Matrix{scalar_type_to_polymake[T]}(undef, 0, size(vr, 2)) : L

    # Rays and Points are homogenized and combined and
    points = homogenize(vr, 1)
    # If some vertices are far vertices, give them a leading 0
    if !isnothing(far_vertices)
        points[far_vertices,1] .= 0
    end

    # Lineality is homogenized
    lineality = homogenize(LM, 0)

    PolyhedralComplex{T}(Polymake.fan.PolyhedralComplex{scalar_type_to_polymake[T]}(
        POINTS = points,
        INPUT_LINEALITY = lineality,
        INPUT_CONES = polyhedra,
    ))
end

# TODO: Only works for this specific case; implement generalization using `iter.Acc`
# Fallback like: PolyhedralFan(itr::AbstractVector{Cone{T}}) where T<:scalar_types
# This makes sure that PolyhedralComplex(maximal_polyhedra(PC)) returns an Oscar PolyhedralComplex,
PolyhedralComplex(iter::SubObjectIterator{Polyhedron{T}}) where T<:scalar_types = PolyhedralComplex{T}(iter.Obj)

###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, PC::PolyhedralComplex{T}) where T<:scalar_types
    try
        ad = ambient_dim(PC)
        print(io, "A polyhedral complex in ambient dimension $(ad)")
        T != fmpq && print(io, " with $T type coefficients")
    catch e
        print(io, "A polyhedral complex without ambient dimension")
    end
end
