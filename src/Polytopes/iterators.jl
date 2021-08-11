# When asking for a property of a `Polymake.BigObject`, e.g. VERTICES, the
# polymake kernel entirely calculates this property before we can even partially
# access this property. Thus, these iterators are constructed with meaningful
# input and perform no calculations themselves during iteration. While this
# implies that the (probably time-critical) calculations have to be performed
# before the actual construction of an iterator (which would otherwise usually
# be executed during the first iteration), this speeds up access and allows for
# useful abstractions and definitions.

struct Polyhedron #a real polymake polyhedron
    pm_polytope::Polymake.BigObject
    boundedness::Symbol # Values: :unknown, :bounded, :unbounded
end

struct Cone #a real polymake polyhedron
    pm_cone::Polymake.BigObject
end

struct Polyhedra
end

struct Points{U} <: AbstractVector{U}
    p::Polymake.Vector{U}
end

Base.IndexStyle(::Type{<:Points}) = IndexLinear()

Base.getindex(po::Points, i::Int64) = po.p[i]

function Base.setindex!(po::Points{U}, val::U, i::Int64) where U
    @boundscheck checkbounds(po.p, i)
    po.p[i] = val
    return val
end

Base.firstindex(::Points) = 1
Base.lastindex(iter::Points) = length(iter)
Base.size(po::Points) = size(po.p)

Points(x...) = Points{Polymake.Rational}(x...)

Points{U}(n::Base.Integer) where U = Points{U}(Polymake.Vector{U}(undef, Polymake.Integer(n)))

function Base.similar(X::Points, ::Type{S}, dims::Dims{1}) where S <: Union{Polymake.Rational, Polymake.Integer}
    return Points{S}(dims...)
end

Base.BroadcastStyle(::Type{<:Points}) = Broadcast.ArrayStyle{Points}()

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{Points}}, ::Type{ElType}) where ElType
    return Points{Polymake.promote_to_pm_type(Vector, ElType)}(axes(bc)...)
end

struct Ray{U} <: AbstractVector{U}
    p::Polymake.Vector{U}
end

Base.IndexStyle(::Type{<:Ray}) = IndexLinear()

Base.getindex(po::Ray, i::Int64) = po.p[i]

function Base.setindex!(po::Ray{U}, val::U, i::Int64) where U
    @boundscheck checkbounds(po.p, i)
    po.p[i] = val
    return val
end

Base.firstindex(::Ray) = 1
Base.lastindex(iter::Ray) = length(iter)
Base.size(po::Ray) = size(po.p)

Ray(x...) = Ray{Polymake.Rational}(x...)

Ray{U}(n::Base.Integer) where U = Ray{U}(Polymake.Vector{U}(undef, Polymake.Integer(n)))

function Base.similar(X::Ray, ::Type{S}, dims::Dims{1}) where S <: Union{Polymake.Rational, Polymake.Integer}
    return Ray{S}(dims...)
end

Base.BroadcastStyle(::Type{<:Ray}) = Broadcast.ArrayStyle{Ray}()

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{Ray}}, ::Type{ElType}) where ElType
    return Ray{Polymake.promote_to_pm_type(Vector, ElType)}(axes(bc)...)
end

struct Rays
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

struct Cones
end

struct PolyhedronOrConeIterator{T} <: AbstractArray{T, 1}
    vertices::Polymake.Matrix{Polymake.Rational}
    faces::Polymake.Array{Polymake.Set{Int64}}
    lineality::Polymake.Matrix{Polymake.Rational}
end

function Base.getindex(iter::PolyhedronOrConeIterator{Polyhedron}, i::Int64)
    @boundscheck checkbounds(iter.faces, i)
    return Polyhedron(Polymake.polytope.Polytope(VERTICES=iter.vertices[[f+1 for f in iter.faces[i]],:],LINEALITY_SPACE = iter.lineality))
end

function Base.getindex(iter::PolyhedronOrConeIterator{Cone}, i::Int64)
    @boundscheck checkbounds(iter.faces, i)
    return Cone(Polymake.polytope.Cone(RAYS=iter.vertices[[f+1 for f in iter.faces[i]],:],LINEALITY_SPACE = iter.lineality))
end

function Base.setindex!(iter::PolyhedronOrConeIterator, val::Polymake.Set{Int64}, i::Int64)
    @boundscheck checkbounds(iter.faces, i)
    iter.faces[i] = val
end

# function setvertices!(iter::PolyhedronOrConeIterator{Polyhedron}, val::Polymake.Matrix{Polymake.Rational}; non_redundant = false)
#     if non_redundant
#         iter.vertices = val
#     else
#         iter.vertices = remove_redundant_rows(val)
#     end
#     return val
# end
#
# function addvertices!(iter::PolyhedronOrConeIterator{Polyhedron}, val::Polymake.Matrix{Polymake.Rational}; non_redundant = false)
#     if non_redundant
#         iter.vertices = [iter.vertices; val]
#     else
#         iter.vertices = remove_redundant_rows([iter.vertices; val])
#     end
#     return val
# end

Base.firstindex(::PolyhedronOrConeIterator) = 1
Base.lastindex(iter::PolyhedronOrConeIterator) = length(iter)
Base.size(iter::PolyhedronOrConeIterator) = (length(iter.faces),)

function Base.show(io::IO, I::PolyhedronOrConeIterator{T}) where T
    print(io, "A collection of $T objects as `$T`")
end

#####################################

struct HalfSpaceIterator{T} <: AbstractArray{T, 1}
    A::Polymake.Matrix{Polymake.Rational}
    b::Polymake.Vector{Polymake.Rational}
end

Base.IndexStyle(::Type{<:HalfSpaceIterator}) = IndexLinear()

function Base.getindex(iter::HalfSpaceIterator{T}, i::Int64) where T
    @boundscheck checkbounds(iter.b, i)
    return T(reshape(iter.A[i, :], 1, :), iter.b[i])
end

function Base.setindex!(iter::HalfSpaceIterator, val::Halfspace, i::Int64)
    @boundscheck checkbounds(iter.b, i)
    iter.A[i, :] = val.a
    iter.b[i] = val.b
    return val
end

Base.firstindex(::HalfSpaceIterator) = 1
Base.lastindex(iter::HalfSpaceIterator) = length(iter)
Base.size(iter::HalfSpaceIterator) = (length(iter.b),)

halfspace_matrix_pair(iter::HalfSpaceIterator) = (A = iter.A, b = iter.b)

function point_matrix(iter::HalfSpaceIterator)
    if !iszero(iter.b)
        throw(ArgumentError("Half-spaces have to be affine in order for point_matrix to be definded"))
    end
    return iter.A
end

HalfSpaceIterator(x...) = HalfSpaceIterator{Halfspace}(x...)

###############################

struct PointIterator{T, U} <: AbstractArray{T, 1}
    m::Polymake.Matrix{U}
end

Base.IndexStyle(::Type{<:PointIterator}) = IndexLinear()

function Base.getindex(iter::PointIterator{T, U}, i::Int64) where {T, U}
    @boundscheck checkbounds(iter.m, i, 1)
    return T{U}(iter.m[i, :])
end

function Base.setindex!(iter::PointIterator{T}, val::T, i::Int64) where T
    @boundscheck checkbounds(iter.m, i, 1)
    iter.m[i, :] = val.p
    return val
end

Base.firstindex(::PointIterator) = 1
Base.lastindex(iter::PointIterator) = length(iter)
Base.size(iter::PointIterator) = (size(iter.m, 1),)

point_matrix(iter::PointIterator) = iter.m

PointIterator{T, U}(vertices::Union{Oscar.fmpz_mat,AbstractMatrix{Oscar.fmpz}, Oscar.fmpq_mat,AbstractMatrix{Oscar.fmpq}}) where {T, U} = PointIterator{T, U}(matrix_for_polymake(vertices))
PointIterator{T, U}(vertices::Array{V, 1}) where {T, U, V<:Union{AbstractVector, T}} = PointIterator{T, U}(cat((v' for v in vertices)...; dims = 1))

PointIterator{T}(x...) where T = PointIterator{T, Polymake.Rational}(x...)
PointIterator(x...) = PointIterator{Points}(x...)

####################

AsTypeIdentities(as::Type{T}) where T<:Union{Points, Polymake.Vector} = Points
AsTypeIdentities(as::Type{T}) where T<:Union{Rays, Ray} = Ray
AsTypeIdentities(as::Type{T}) where T<:Union{Halfspace, Halfspaces} = Halfspace
AsTypeIdentities(as::Type{T}) where T<:Union{Polyhedron, Polyhedra} = Polyhedron
AsTypeIdentities(as::Type{T}) where T<:Pair = Pair{Polymake.Matrix, Polymake.Rational}
AsTypeIdentities(as::Type{T}) where T<:Union{Cone, Cones} = Cone

####################

matrix_for_polymake(iter::PointIterator) = iter.m
