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

function Base.getindex(iter::PolyhedronFaceIterator, i::Int64)
    @boundscheck checkbounds(iter.faces, i)
    return iterate(iter, i)[1]
end

function Base.show(io::IO, I::PolyhedronFaceIterator{T}) where T
    print(io, "A collection of faces as `$T`")
end

#####################################

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

function Base.getindex(iter::PolyhedronFacetIterator, i::Int64)
    @boundscheck checkbounds(iter.b, i)
    return iterate(iter, i)[1]
end

halfspace_matrix_pair(iter::PolyhedronFacetIterator) = (A = iter.A, b = iter.b)

PolyhedronFacetIterator{T}(x...) where T<:Halfspaces = PolyhedronFacetIterator{Halfspace}(x...)
PolyhedronFacetIterator{T}(x...) where T<:Polyhedra = PolyhedronFacetIterator{Polyhedron}(x...)

function Base.show(io::IO, I::PolyhedronFacetIterator{T}) where T
    print(io, "A collection of facets as `$T`")
end

###############################

struct PointIterator{T, U}
    m::Polymake.Matrix{U}
end

function Base.iterate(iter::PointIterator, index::Int64 = 1)
    if (index <= size(iter.m, 1))
        return (iter.m[index, :], index +1)
    end
    return nothing
end

Base.eltype(::Type{PointIterator{T, U}}) where {T, U} = Polymake.Vector{U}
Base.length(iter::PointIterator) = size(iter.m, 1)

function Base.getindex(iter::PointIterator, i::Int64)
    @boundscheck checkbounds(iter.m, i, 1)
    return iterate(iter, i)[1]
end

point_matrix(iter::PointIterator) = iter.m

function Base.show(io::IO, I::PointIterator{T, U}) where {T, U}
    print(io, "A collection of faces as `$T{$U}`")
end

####################

AsTypeIdentitiesP(as::Type{T}) where T<:Points = Polymake.Vector
AsTypeIdentitiesF(as::Type{T}) where T<:Halfspaces = Halfspace
AsTypeIdentitiesF(as::Type{T}) where T<:Polyhedra = Polyhedron
AsTypeIdentitiesFD(as::Type{T}) where T<:Polyhedra = Polyhedron
