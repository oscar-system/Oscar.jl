# struct PointIterator
#     o::Polymake.BigObject
#     p::Symbol
# end
#
# # VERTICES
# function _iterate(o::Polymake.BigObject, ::Val{:VERTICES}, index)
#     vertices = o.VERTICES
#     while index <= size(vertices, 1)
#         if iszero(vertices[index, 1])
#             index += 1
#         else
#             return (vertices[index, 2:end], index + 1)
#         end
#     end
#     return nothing
# end
#
# _length(o::Polymake.BigObject, ::Val{:VERTICES}) = o.N_VERTICES - length(o.FAR_FACE)
#
# # LATTICE_POINTS
# function _iterate(o::Polymake.BigObject, ::Val{:LATTICE_POINTS_GENERATORS}, index)
#     lp = o.LATTICE_POINTS_GENERATORS[1]
#     while index <= size(lp, 1)
#         if iszero(lp[index, 1])
#             index += 1
#         else
#             return (lp[index, 2:end], index + 1)
#         end
#     end
#     return nothing
# end
#
# _length(o::Polymake.BigObject, ::Val{:LATTICE_POINTS_GENERATORS}) = size(o.LATTICE_POINTS_GENERATORS[1],2)
#
# function Base.iterate(iter::PointIterator, index::Int64 = 1)
#     return _iterate(iter.o, Val(iter.p), index)
# end
#
# Base.eltype(::Type{PointIterator}) = Polymake.Vector{Polymake.Rational}
# Base.length(iter::PointIterator) = _length(iter.o, Val(iter.p))

struct Polyhedra
end
struct Points
end
@doc Markdown.doc"""

    Halfspaces

Dummy type used for specifying the desired output format.
One halfspace `H(a,b)` is given by a vector `a` and a value `b` such that
$$H(a,b) = \{ x | ax ≤ b \}.$$
"""
struct Halfspaces
end

struct Halfspace
    a::Polymake.Matrix{Polymake.Rational}
    b::Polymake.Rational
end

struct PolyhedronFaceIterator{T}
    p::Polymake.BigObject
    faces::Polymake.Array{Polymake.Set{Int64}}
end

function Base.iterate(iter::PolyhedronFaceIterator, index::Int64 = 1)
    if index > length(iter.faces)
        return nothing
    end
    p = Polyhedron(Polymake.polytope.Polytope(VERTICES=iter.p.VERTICES[[f+1 for f in iter.faces[index]],:],LINEALITY_SPACE = iter.p.LINEALITY_SPACE))
    return (p,index+1)
end

Base.eltype(::Type{PolyhedronFaceIterator{T}}) where T = T
Base.length(iter::PolyhedronFaceIterator) = length(iter.faces)

#####################################

# struct PolyhedronFacetIterator{T}
#     facets::Polymake.Matrix{Polymake.Rational}
# end

struct PolyhedronFacetIterator{T}
    A::AbstractMatrix{Polymake.Rational}
    b::AbstractVector{Polymake.Rational}
end

function Base.iterate(iter::PolyhedronFacetIterator{T}, index::Int64 = 1) where T
    if (length(iter.b) < index)
        return nothing
    end
    return (T(reshape(iter.A[index, :], 1, :), iter.b[index]), index + 1)
end

Base.eltype(::Type{PolyhedronFacetIterator{T}}) where T = T
Base.length(iter::PolyhedronFacetIterator) = length(iter.b)

halfspace_matrix_pair(iter::PolyhedronFacetIterator) = (A = iter.A, b = iter.b)

###############################

struct PointIterator
    m::Polymake.Matrix{Polymake.Rational}
end

function Base.iterate(iter::PointIterator, index::Int64 = 1)
    if (index <= size(iter.m, 1))
        return (iter.m[index, :], index +1)
    end
    return nothing
end

Base.eltype(::Type{PointIterator}) = Polymake.Vector{Polymake.Rational}
Base.length(iter::PointIterator) = size(iter.m, 1)

function Base.getindex(iter::PointIterator, i::Int64)
    @boundscheck checkbounds(iter.m, i, 1)
    return iterate(iter)[1]
end

point_matrix(iter::PointIterator) = iter.m
