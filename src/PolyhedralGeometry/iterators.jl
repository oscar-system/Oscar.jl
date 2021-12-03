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

abstract type Halfspace end

@doc Markdown.doc"""
    Halfspace(a, b)

One halfspace `H(a,b)` is given by a vector `a` and a value `b` such that
$$H(a,b) = \{ x | ax ≤ b \}.$$
"""
struct AffineHalfspace <: Halfspace
    a::Polymake.Vector{Polymake.Rational}
    b::Polymake.Rational
end

AffineHalfspace(a::Union{MatElem, AbstractMatrix}, b=0) = AffineHalfspace(vec(a), b)

AffineHalfspace(a) = AffineHalfspace(a, 0)

Halfspace(a, b) = AffineHalfspace(a, b)

negbias(H::AffineHalfspace) = H.b

struct LinearHalfspace <: Halfspace
    a::Polymake.Vector{Polymake.Rational}
end

LinearHalfspace(a::Union{MatElem, AbstractMatrix}) = LinearHalfspace(vec(a))

Halfspace(a) = LinearHalfspace(a)

negbias(H::LinearHalfspace) = 0

abstract type Hyperplane end

@doc Markdown.doc"""
    AffineHyperplane(a, b)

One hyperplane `H(a,b)` is given by a vector `a` and a value `b` such that
$$H(a,b) = \{ x | ax = b \}.$$
"""
struct AffineHyperplane <: Hyperplane
    a::Polymake.Vector{Polymake.Rational}
    b::Polymake.Rational
end

AffineHyperplane(a::Union{MatElem, AbstractMatrix}, b) = AffineHyperplane(vec(a), b)

AffineHyperplane(a) = AffineHyperplane(a, 0)

Hyperplane(a, b) = AffineHyperplane(a, b)

negbias(H::AffineHyperplane) = H.b

struct LinearHyperplane <: Hyperplane
    a::Polymake.Vector{Polymake.Rational}
end

LinearHyperplane(a::Union{MatElem, AbstractMatrix}) = LinearHyperplane(vec(a))

Hyperplane(a) = LinearHyperplane(a)

negbias(H::LinearHyperplane) = 0

# TODO: abstract notion of equality
Base.:(==)(x::AffineHalfspace, y::AffineHalfspace) = x.a == y.a && x.b == y.b

Base.:(==)(x::LinearHalfspace, y::LinearHalfspace) = x.a == y.a

Base.:(==)(x::AffineHyperplane, y::AffineHyperplane) = x.a == y.a && x.b == y.b

Base.:(==)(x::LinearHyperplane, y::LinearHyperplane) = x.a == y.a

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

ray_incidences(iter::SubObjectIterator{<:Union{Cone, Polyhedron}}) = _ray_incidences(Val(iter.Acc), iter.Obj; iter.options...)
_ray_incidences(::Any, ::Polymake.BigObject) = throw(ArgumentError("Incidence Matrix resp. rays not defined in this context."))

vertex_incidences(iter::SubObjectIterator{Polyhedron}) = _vertex_incidences(Val(iter.Acc), iter.Obj; iter.options...)
_vertex_incidences(::Any, ::Polymake.BigObject) = throw(ArgumentError("Incidence Matrix resp. vertices not defined in this context."))

linear_inequality_matrix(iter::SubObjectIterator{<:Union{Cone, Polyhedron, Halfspace}}) = matrix(QQ, Matrix{fmpq}(_linear_inequality_matrix(Val(iter.Acc), iter.Obj; iter.options...)))
_linear_inequality_matrix(::Any, ::Polymake.BigObject) = throw(ArgumentError("Linear Inequality Matrix not defined in this context."))

affine_inequality_matrix(iter::SubObjectIterator{<:Union{Halfspace, Hyperplane, Pair, Polyhedron}}) = matrix(QQ, Matrix{fmpq}(_affine_inequality_matrix(Val(iter.Acc), iter.Obj; iter.options...)))
_affine_inequality_matrix(::Any, ::Polymake.BigObject) = throw(ArgumentError("Affine Inequality Matrix not defined in this context."))

linear_equation_matrix(iter::SubObjectIterator{<:Union{Halfspace, Hyperplane, Pair, Polyhedron}}) = matrix(QQ, Matrix{fmpq}(_linear_equation_matrix(Val(iter.Acc), iter.Obj; iter.options...)))
_linear_equation_matrix(::Any, ::Polymake.BigObject) = throw(ArgumentError("Linear Equation Matrix not defined in this context."))

affine_equation_matrix(iter::SubObjectIterator{<:Union{Halfspace, Hyperplane, Pair, Polyhedron}}) = matrix(QQ, Matrix{fmpq}(_affine_equation_matrix(Val(iter.Acc), iter.Obj; iter.options...)))
_affine_equation_matrix(::Any, ::Polymake.BigObject) = throw(ArgumentError("Affine Equation Matrix not defined in this context."))

function matrix_for_polymake(iter::SubObjectIterator)
    if hasmethod(_matrix_for_polymake, Tuple{Val{iter.Acc}})
        return _matrix_for_polymake(Val(iter.Acc))(Val(iter.Acc), iter.Obj; iter.options...)
    else
        throw(ArgumentError("Matrix for Polymake not defined in this context."))
    end
end

function linear_matrix_for_polymake(iter::SubObjectIterator{<:Union{Halfspace, Hyperplane, Pair, Polyhedron}})
    if hasmethod(_linear_matrix_for_polymake, Tuple{Val{iter.Acc}})
        return _linear_matrix_for_polymake(Val(iter.Acc))(Val(iter.Acc), iter.Obj; iter.options...)
    elseif hasmethod(_affine_matrix_for_polymake, Tuple{Val{iter.Acc}})
        res = _affine_matrix_for_polymake(Val(iter.Acc))(Val(iter.Acc), iter.Obj; iter.options...)
        iszero(res[:, 1]) || throw(ArgumentError("Input not linear."))
        return res[:, 2:end]
    end
    throw(ArgumentError("Linear Matrix for Polymake not defined in this context."))
end
_linear_matrix_for_polymake(::Any, ::Polymake.BigObject) = throw(ArgumentError("Linear Matrix for Polymake not defined in this context."))

function affine_matrix_for_polymake(iter::SubObjectIterator{<:Union{Halfspace, Hyperplane, Pair, Polyhedron}})
    if hasmethod(_affine_matrix_for_polymake, Tuple{Val{iter.Acc}})
        return _affine_matrix_for_polymake(Val(iter.Acc))(Val(iter.Acc), iter.Obj; iter.options...)
    elseif hasmethod(_linear_matrix_for_polymake, Tuple{Val{iter.Acc}})
        return homogenize(_linear_matrix_for_polymake(Val(iter.Acc))(Val(iter.Acc), iter.Obj; iter.options...), 0)
    end
    throw(ArgumentError("Affine Matrix for Polymake not defined in this context."))
end
_affine_matrix_for_polymake(::Any, ::Polymake.BigObject) = throw(ArgumentError("Affine Matrix for Polymake not defined in this context."))

function halfspace_matrix_pair(iter::SubObjectIterator{<:Union{Halfspace, Hyperplane, Pair, Polyhedron}})
    h = affine_matrix_for_polymake(iter)
    return (A = matrix(QQ, Matrix{fmpq}(h[:, 2:end])), b = -h[:, 1])
end
_halfspace_matrix_pair(::Any, ::Polymake.BigObject) = throw(ArgumentError("Halfspace Matrix Pair not defined in this context."))

point_matrix(iter::SubObjectIterator{<:AbstractVector{Polymake.Rational}}) = matrix(QQ, Matrix{fmpq}(_point_matrix(Val(iter.Acc), iter.Obj; iter.options...)))
point_matrix(iter::SubObjectIterator{<:AbstractVector{Polymake.Integer}}) = matrix(ZZ, _point_matrix(Val(iter.Acc), iter.Obj; iter.options...))
_point_matrix(::Any, ::Polymake.BigObject) = throw(ArgumentError("Point Matrix not defined in this context."))

vector_matrix(iter::SubObjectIterator{<:AbstractVector{Polymake.Rational}}) = matrix(QQ, Matrix{fmpq}(_vector_matrix(Val(iter.Acc), iter.Obj; iter.options...)))
vector_matrix(iter::SubObjectIterator{<:AbstractVector{Polymake.Integer}}) = matrix(ZZ, _vector_matrix(Val(iter.Acc), iter.Obj; iter.options...))
_vector_matrix(::Any, ::Polymake.BigObject) = throw(ArgumentError("Vector Matrix not defined in this context."))

generator_matrix(iter::SubObjectIterator{<:AbstractVector{Polymake.Rational}}) = matrix(QQ, Matrix{fmpq}(_generator_matrix(Val(iter.Acc), iter.Obj; iter.options...)))
generator_matrix(iter::SubObjectIterator{<:AbstractVector{Polymake.Integer}}) = matrix(ZZ, _generator_matrix(Val(iter.Acc), iter.Obj; iter.options...))
_generator_matrix(::Any, ::Polymake.BigObject) = throw(ArgumentError("Generator Matrix not defined in this context."))
