###############################################################################
###############################################################################
### Definition and constructors
###############################################################################
###############################################################################

struct PolyhedralComplex{T} <: PolyhedralObject{T}
    pm_complex::Polymake.BigObject
    parent_field::Field
    
    PolyhedralComplex{T}(pm::Polymake.BigObject, p::Field) where T<:scalar_types = new{T}(pm, p)
end


polyhedral_complex(p::Polymake.BigObject) = PolyhedralComplex{detect_scalar_type(PolyhedralComplex, p)}(p)

pm_object(pc::PolyhedralComplex) = pc.pm_complex


@doc raw"""
    polyhedral_complex(::T, polyhedra, vr, far_vertices, L) where T<:scalar_types

# Arguments
- `T`: Type of scalar to use, defaults to `QQFieldElem`.
- `polyhedra::IncidenceMatrix`: An incidence matrix; there is a 1 at position
  (i,j) if the ith polytope contains point j and 0 otherwise.
- `vr::AbstractCollection[PointVector]`: The points whose convex hulls make up
  the polyhedral complex. This matrix also contains the far vertices.
- `far_vertices::Vector{Int}`: Vector containing the indices of the rows
  corresponding to the far vertices in `vr`.
- `L::AbstractCollection[RayVector]`: Generators of the lineality space of the
  polyhedral complex.

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

julia> PC = polyhedral_complex(IM, vr)
Polyhedral complex in ambient dimension 2
```

Polyhedral complex with rays and lineality:
```jldoctest
julia> VR = [0 0 0; 1 0 0; 0 1 0; -1 0 0];

julia> IM = IncidenceMatrix([[1,2,3],[1,3,4]]);

julia> far_vertices = [2,3,4];

julia> L = [0 0 1];

julia> PC = polyhedral_complex(IM, VR, far_vertices, L)
Polyhedral complex in ambient dimension 3

julia> lineality_dim(PC)
1
```
"""
function polyhedral_complex(f::Union{Type{T}, Field},
                            polyhedra::IncidenceMatrix, 
                            vr::AbstractCollection[PointVector], 
                            far_vertices::Union{Vector{Int}, Nothing} = nothing, 
                            L::Union{AbstractCollection[RayVector], Nothing} = nothing;
                            non_redundant::Bool = false
                            ) where T<:scalar_types
    parent_field, scalar_type = _determine_parent_and_scalar(f, vr, L)
    LM = isnothing(L) || isempty(L) ? Polymake.Matrix{_scalar_type_to_polymake(scalar_type)}(undef, 0, size(vr, 2)) : L

    # Rays and Points are homogenized and combined and
    points = homogenize(vr, 1)
    # If some vertices are far vertices, give them a leading 0
    if !isnothing(far_vertices)
        points[far_vertices,1] .= 0
    end

    # Lineality is homogenized
    lineality = homogenize(LM, 0)

    if non_redundant
        return PolyhedralComplex{scalar_type}(Polymake.fan.PolyhedralComplex{_scalar_type_to_polymake(scalar_type)}(
            VERTICES = points,
            LINEALITY_SPACE = lineality,
            MAXIMAL_CONES = polyhedra,
        ), parent_field)
    else
        return PolyhedralComplex{scalar_type}(Polymake.fan.PolyhedralComplex{_scalar_type_to_polymake(scalar_type)}(
            POINTS = points,
            INPUT_LINEALITY = lineality,
            INPUT_CONES = polyhedra,
        ), parent_field)
    end
end
# default scalar type: `QQFieldElem`
polyhedral_complex(polyhedra::IncidenceMatrix, 
                vr::AbstractCollection[PointVector], 
                far_vertices::Union{Vector{Int}, Nothing} = nothing, 
                L::Union{AbstractCollection[RayVector], Nothing} = nothing;
                non_redundant::Bool = false) =
  polyhedral_complex(QQFieldElem, polyhedra, vr, far_vertices, L; non_redundant=non_redundant)

# TODO: Only works for this specific case; implement generalization using `iter.Acc`
# Fallback like: PolyhedralFan(itr::AbstractVector{Cone{T}}) where T<:scalar_types
# This makes sure that polyhedral_complex(maximal_polyhedra(PC)) returns an Oscar PolyhedralComplex,
polyhedral_complex(iter::SubObjectIterator{Polyhedron{T}}) where T<:scalar_types = PolyhedralComplex{T}(iter.Obj)

###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, PC::PolyhedralComplex{T}) where T<:scalar_types
    try
        ad = ambient_dim(PC)
        print(io, "Polyhedral complex in ambient dimension $(ad)")
        T != QQFieldElem && print(io, " with $T type coefficients")
    catch e
        print(io, "Polyhedral complex without ambient dimension")
    end
end
