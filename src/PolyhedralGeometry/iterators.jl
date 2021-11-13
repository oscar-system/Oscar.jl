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

matrix_for_polymake(iter::VectorIterator) = iter.m

####################

struct SubObjectIterator{T} <: AbstractVector{T}
    Obj::Polymake.BigObject
    Acc::Function
    n::Base.Integer
    options::NamedTuple
end

SubObjectIterator{T}(Obj::Polymake.BigObject, Acc::Function, n::Base.Integer) where T = SubObjectIterator{T}(Obj, Acc, n, NamedTuple())

Base.IndexStyle(::Type{<:SubObjectIterator}) = IndexLinear()

function Base.getindex(iter::SubObjectIterator{T}, i::Base.Integer) where T
    @boundscheck 1 <= i && i <= iter.n
    return iter.Acc(T, iter.Obj, i; iter.options...)
end

Base.firstindex(::SubObjectIterator) = 1
Base.lastindex(iter::SubObjectIterator) = length(iter)
Base.size(iter::SubObjectIterator) = (iter.n,)

ray_incidences(iter::SubObjectIterator) = _ray_incidences(Val(iter.Acc), iter.Obj; iter.options...)
_ray_incidences(::Any, ::Polymake.BigObject) = throw(ArgumentError("Incidence Matrix resp. rays not defined in this context."))

vertex_incidences(iter::SubObjectIterator) = _vertex_incidences(Val(iter.Acc), iter.Obj; iter.options...)
_vertex_incidences(::Any, ::Polymake.BigObject) = throw(ArgumentError("Incidence Matrix resp. vertices not defined in this context."))

inequality_matrix(iter::SubObjectIterator) = matrix(QQ, Matrix{fmpq}(_inequality_matrix(Val(iter.Acc), iter.Obj; iter.options...)))
_inequality_matrix(::Any, ::Polymake.BigObject) = throw(ArgumentError("Inequality Matrix not defined in this context."))

affine_inequality_matrix(iter::SubObjectIterator) = matrix(QQ, Matrix{fmpq}(_affine_inequality_matrix(Val(iter.Acc), iter.Obj; iter.options...)))
_affine_inequality_matrix(::Any, ::Polymake.BigObject) = throw(ArgumentError("Inequality Matrix not defined in this context."))

equation_matrix(iter::SubObjectIterator) = matrix(QQ, Matrix{fmpq}(_equation_matrix(Val(iter.Acc), iter.Obj; iter.options...)))
_equation_matrix(::Any, ::Polymake.BigObject) = throw(ArgumentError("Equation Matrix not defined in this context."))

affine_equation_matrix(iter::SubObjectIterator) = matrix(QQ, Matrix{fmpq}(_affine_equation_matrix(Val(iter.Acc), iter.Obj; iter.options...)))
_affine_equation_matrix(::Any, ::Polymake.BigObject) = throw(ArgumentError("Affine Equation Matrix not defined in this context."))

matrix_for_polymake(iter::SubObjectIterator) = _matrix_for_polymake(Val(iter.Acc), iter.Obj; iter.options...)
_matrix_for_polymake(::Any, ::Polymake.BigObject) = throw(ArgumentError("Matrix for Polymake not defined in this context."))

affine_matrix_for_polymake(iter::SubObjectIterator) = _affine_matrix_for_polymake(Val(iter.Acc), iter.Obj; iter.options...)
_affine_matrix_for_polymake(::Any, ::Polymake.BigObject) = throw(ArgumentError("Affine Matrix for Polymake not defined in this context."))

function halfspace_matrix_pair(iter::SubObjectIterator)
    h = _halfspace_matrix_pair(Val(iter.Acc), iter.Obj; iter.options...)
    return (A = matrix(QQ, Matrix{fmpq}(h[1])), b = h[2])
end
_halfspace_matrix_pair(::Any, ::Polymake.BigObject) = throw(ArgumentError("Halfspace Matrix Pair not defined in this context."))

vector_matrix(iter::SubObjectIterator) = matrix(QQ, Matrix{fmpq}(_vector_matrix(Val(iter.Acc), iter.Obj; iter.options...)))
_vector_matrix(::Any, ::Polymake.BigObject) = throw(ArgumentError("Vector Matrix not defined in this context."))

generator_matrix(iter::SubObjectIterator) = matrix(QQ, Matrix{fmpq}(_generator_matrix(Val(iter.Acc), iter.Obj; iter.options...)))
_generator_matrix(::Any, ::Polymake.BigObject) = throw(ArgumentError("Generator Matrix not defined in this context."))
