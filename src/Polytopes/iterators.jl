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

struct PointVector{U} <: AbstractVector{U}
    p::Polymake.Vector
    PointVector{U}(p) where U = new(Polymake.Vector{Polymake.to_cxx_type(U)}(p))
end

Base.IndexStyle(::Type{<:PointVector}) = IndexLinear()

Base.getindex(po::PointVector, i::Base.Integer) = po.p[i]

function Base.setindex!(po::PointVector{U}, val, i::Base.Integer) where U
    @boundscheck checkbounds(po.p, i)
    po.p[i] = val
    return val
end

Base.firstindex(::PointVector) = 1
Base.lastindex(iter::PointVector) = length(iter)
Base.size(po::PointVector) = size(po.p)

PointVector(x...) = PointVector{Polymake.Rational}(x...)

PointVector{U}(n::Base.Integer) where U = PointVector{U}(Polymake.Vector{U}(undef, Polymake.Integer(n)))

function Base.similar(X::PointVector, ::Type{S}, dims::Dims{1}) where S <: Union{Polymake.Rational, Polymake.Integer, Int64, Float64}
    return PointVector{S}(dims...)
end

Base.BroadcastStyle(::Type{<:PointVector}) = Broadcast.ArrayStyle{PointVector}()

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{PointVector}}, ::Type{ElType}) where ElType
    return PointVector{Polymake.promote_to_pm_type(Vector, ElType)}(axes(bc)...)
end

struct RayVector{U} <: AbstractVector{U}
    p::Polymake.Vector
    RayVector{U}(p) where U = new(Polymake.Vector{Polymake.to_cxx_type(U)}(p))
end

Base.IndexStyle(::Type{<:RayVector}) = IndexLinear()

Base.getindex(po::RayVector, i::Base.Integer) = po.p[i]

function Base.setindex!(po::RayVector{U}, val, i::Base.Integer) where U
    @boundscheck checkbounds(po.p, i)
    po.p[i] = val
    return val
end

Base.firstindex(::RayVector) = 1
Base.lastindex(iter::RayVector) = length(iter)
Base.size(po::RayVector) = size(po.p)

RayVector(x...) = RayVector{Polymake.Rational}(x...)

RayVector{U}(n::Base.Integer) where U = RayVector{U}(Polymake.Vector{U}(undef, Polymake.Integer(n)))

function Base.similar(X::RayVector, ::Type{S}, dims::Dims{1}) where S <: Union{Polymake.Rational, Polymake.Integer, Int64, Float64}
    return RayVector{S}(dims...)
end

Base.BroadcastStyle(::Type{<:RayVector}) = Broadcast.ArrayStyle{RayVector}()

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{RayVector}}, ::Type{ElType}) where ElType
    return RayVector{Polymake.promote_to_pm_type(Vector, ElType)}(axes(bc)...)
end

@doc Markdown.doc"""
    Halfspace(a, b)

One halfspace `H(a,b)` is given by a vector `a` and a value `b` such that
$$H(a,b) = \{ x | ax ≤ b \}.$$
"""
struct Halfspace
    a::Polymake.Vector{Polymake.Rational}
    b::Polymake.Rational
end

Halfspace(a::Union{MatElem, AbstractMatrix}, b=0) = Halfspace(vec(a), b)

Halfspace(a) = Halfspace(a, 0)

@doc Markdown.doc"""
    Hyperplane(a, b)

One hyperplane `H(a,b)` is given by a vector `a` and a value `b` such that
$$H(a,b) = \{ x | ax = b \}.$$
"""
struct Hyperplane
    a::Polymake.Vector{Polymake.Rational}
    b::Polymake.Rational
end

Hyperplane(a::Union{MatElem, AbstractMatrix}, b) = Hyperplane(vec(a), b)

Hyperplane(a) = Hyperplane(a, 0)

# TODO: abstract notion of equality
Base.:(==)(x::Halfspace, y::Halfspace) = x.a == y.a && x.b == y.b

Base.:(==)(x::Hyperplane, y::Hyperplane) = x.a == y.a && x.b == y.b

struct PolyhedronOrConeIterator{T} <: AbstractVector{T}
    vertices::Polymake.Matrix{Polymake.Rational}
    faces::Polymake.Array{Polymake.Set{Polymake.to_cxx_type(Int64)}}
    lineality::Polymake.Matrix{Polymake.Rational}
end

function Base.getindex(iter::PolyhedronOrConeIterator{Polyhedron}, i::Base.Integer)
    @boundscheck checkbounds(iter.faces, i)
    return Polyhedron(Polymake.polytope.Polytope(VERTICES=iter.vertices[[f for f in iter.faces[i]],:],LINEALITY_SPACE = iter.lineality))
end

function Base.getindex(iter::PolyhedronOrConeIterator{Cone}, i::Base.Integer)
    @boundscheck checkbounds(iter.faces, i)
    return Cone(Polymake.polytope.Cone(RAYS=iter.vertices[[f for f in iter.faces[i]],:],LINEALITY_SPACE = iter.lineality))
end

function Base.setindex!(iter::PolyhedronOrConeIterator, val::Polymake.Set{Polymake.to_cxx_type(Int64)}, i::Base.Integer)
    @boundscheck checkbounds(iter.faces, i)
    iter.faces[i, :] = Polymake.spzeros(size(iter.faces, 2))
    for j in val
        iter.faces[i, j] = true
    end
    return val
end

Base.firstindex(::PolyhedronOrConeIterator) = 1
Base.lastindex(iter::PolyhedronOrConeIterator) = length(iter)
Base.size(iter::PolyhedronOrConeIterator) = (length(iter.faces),)

# TODO: function incidence_matrix(::PolyhedronOrConeIterator) after merge of combine-incidencematrix

function Base.show(io::IO, I::PolyhedronOrConeIterator{T}) where T
    print(io, "A collection of `$T` objects")
end

#####################################

struct HalfspaceIterator{T} <: AbstractVector{T}
    A::Polymake.Matrix{Polymake.Rational}
    b::Polymake.Vector{Polymake.Rational}
end

Base.IndexStyle(::Type{<:HalfspaceIterator}) = IndexLinear()

function Base.getindex(iter::HalfspaceIterator{T}, i::Base.Integer) where T
    @boundscheck checkbounds(iter.b, i)
    return T(reshape(iter.A[i, :], 1, :), iter.b[i])
end

function Base.getindex(iter::HalfspaceIterator{Cone}, i::Base.Integer)
    @boundscheck checkbounds(iter.b, i)
    return cone_from_inequalities(reshape(iter.A[i, :], 1, :))
end

function Base.setindex!(iter::HalfspaceIterator, val::Halfspace, i::Base.Integer)
    @boundscheck checkbounds(iter.b, i)
    iter.A[i, :] = val.a
    iter.b[i] = val.b
    return val
end

Base.firstindex(::HalfspaceIterator) = 1
Base.lastindex(iter::HalfspaceIterator) = length(iter)
Base.size(iter::HalfspaceIterator) = (length(iter.b),)

halfspace_matrix_pair(iter::HalfspaceIterator) = (A = matrix(QQ, Matrix{fmpq}(iter.A)), b = Vector{fmpq}(iter.b))

# Affine Halfspaces

HalfspaceIterator{T}(A::AbstractMatrix) where T = HalfspaceIterator{T}(A, zeros(size(A, 1)))

function matrix(iter::HalfspaceIterator)
    if !iszero(iter.b)
        throw(ArgumentError("Description has to be non-affine in order for `matrix` to be defined"))
    end
    return matrix(QQ, Matrix{fmpq}(iter.A))
end

function matrix_for_polymake(iter::HalfspaceIterator)
    if !iszero(iter.b)
        throw(ArgumentError("Description has to be non-affine in order for `matrix_for_polymake` to be defined"))
    end
    return iter.A
end

affine_matrix_for_polymake(iter::HalfspaceIterator) = hcat(-iter.b, iter.A)

HalfspaceIterator(x...) = HalfspaceIterator{Halfspace}(x...)

###############################

struct VectorIterator{T} <: AbstractVector{T}
    m::AbstractMatrix
    VectorIterator{PointVector{U}}(m) where U = new(Polymake.Matrix{Polymake.to_cxx_type(U)}(m))
    VectorIterator{RayVector{U}}(m) where U = new(Polymake.Matrix{Polymake.to_cxx_type(U)}(m))
end

Base.IndexStyle(::Type{<:VectorIterator}) = IndexLinear()

function Base.getindex(iter::VectorIterator{T}, i::Base.Integer) where T
    @boundscheck checkbounds(iter.m, i, 1)
    return T(iter.m[i, :])
end

function Base.setindex!(iter::VectorIterator{T}, val, i::Base.Integer) where T
    @boundscheck checkbounds(iter.m, i, 1)
    iter.m[i, :] = val
    return val
end

Base.firstindex(::VectorIterator) = 1
Base.lastindex(iter::VectorIterator) = length(iter)
Base.size(iter::VectorIterator) = (size(iter.m, 1),)

matrix(iter::VectorIterator{T}) where {S<:Base.Integer,T<:Union{PointVector{S}, RayVector{S}}} = matrix(ZZ, iter.m)

matrix(iter::VectorIterator{T}) where {S,T<:Union{PointVector{S}, RayVector{S}}} = matrix(QQ, Matrix{fmpq}(iter.m))

VectorIterator{RayVector{S}}(vertices::Union{Oscar.fmpz_mat,AbstractMatrix{Oscar.fmpz}, Oscar.fmpq_mat,AbstractMatrix{Oscar.fmpq}}) where S = VectorIterator{RayVector{S}}(matrix_for_polymake(vertices))
VectorIterator{PointVector{S}}(vertices::Union{Oscar.fmpz_mat,AbstractMatrix{Oscar.fmpz}, Oscar.fmpq_mat,AbstractMatrix{Oscar.fmpq}}) where S = VectorIterator{PointVector{S}}(matrix_for_polymake(vertices))
VectorIterator{T}(vertices::Array{V, 1}) where {T, V<:Union{AbstractVector, T}} = VectorIterator{T}(cat((v' for v in vertices)...; dims = 1))

VectorIterator(x...) = VectorIterator{PointVector{Polymake.Rational}}(x...)

####################

matrix_for_polymake(iter::VectorIterator) = iter.m
