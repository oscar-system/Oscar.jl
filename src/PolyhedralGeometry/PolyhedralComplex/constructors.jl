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
    PolyhedralComplex(Polytopes, Points, Rays, Lineality)

# Arguments
- `Polytopes::IncidenceMatrix`: An incidence matrix; there is a 1 at position
  (i,j) if the ith polytope contains point j and 0 otherwise.
- `Points::Matrix`: The points whose convex hulls make up the polyhedral
  complex.

A polyhedral complex formed from points, rays, and lineality combined into
polyhedra indicated by an incidence matrix, where the columns represent the
points and the rows represent the polyhedra.

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2,3],[1,3,4]])
2×4 IncidenceMatrix
[1, 2, 3]
[1, 3, 4]


julia> V = [0 0; 1 0; 1 1; 0 1]
4×2 Matrix{Int64}:
 0  0
 1  0
 1  1
 0  1

julia> PC = PolyhedralComplex(IM, V)
A polyhedral complex in ambient dimension 2
```
"""
function PolyhedralComplex(Incidence::IncidenceMatrix, V::Union{SubObjectIterator{<:PointVector}, Oscar.MatElem,AbstractMatrix}, R::Union{SubObjectIterator{<:RayVector}, Oscar.MatElem,AbstractMatrix, Nothing} = nothing, L::Union{SubObjectIterator{<:RayVector}, Oscar.MatElem, AbstractMatrix, Nothing} = nothing)
    VM = matrix_for_polymake(V)
    RM = isnothing(R) || isempty(R) ? Polymake.Matrix{Polymake.Rational}(undef, 0, size(VM, 2)) : matrix_for_polymake(R)
    LM = isnothing(L) || isempty(L) ? Polymake.Matrix{Polymake.Rational}(undef, 0, size(VM, 2)) : matrix_for_polymake(L)

    # Rays and Points are homogenized and combined and
    # Lineality is homogenized
    points = stack(homogenize(VM, 1), homogenize(RM, 0))
    lineality = homogenize(LM, 0)

    PolyhedralComplex(Polymake.fan.PolyhedralComplex{Polymake.Rational}(
        POINTS = matrix_for_polymake(points),
        INPUT_LINEALITY = matrix_for_polymake(lineality),
        INPUT_CONES = Incidence,
    ))
end


###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, PC::PolyhedralComplex)
    ad = ambient_dim(PC)
    if ad == -1.0
        print(io, "A polyhedral complex without ambient dimension")
    else
        print(io, "A polyhedral complex in ambient dimension $(ad)")
    end
end

