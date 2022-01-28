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

@doc Markdown.doc"""
    SubObjectIterator(Obj, Acc, n, [options])

An iterator over a designated property of `Obj::Polymake.BigObject`.

`Acc::Function` will be used internally for `getindex`. Further this uniquely
determines the context the iterator operates in, allowing to extend specific
methods like `point_matrix`.

The length of the iterator is hard set with `n::Int`. This is because it is
fixed and to avoid redundant computations: when data has to be pre-processed
before creating a `SubObjectIterator`, the length can usually easily be derived.

Additional data required for specifying the property can be given using
`options::NamedTuple`. A typical example for this is `dim` in the context of
`facets`. The `NamedTuple` is passed to `Acc` (and the specific methods) as
keyword arguments.
"""
struct SubObjectIterator{T} <: AbstractVector{T}
    Obj::Polymake.BigObject
    Acc::Function
    n::Int
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

# Incidence matrices
for (sym, name) in (("ray_indices", "Incidence Matrix resp. rays"), ("vertex_indices", "Incidence Matrix resp. vertices"))
    M = Symbol(sym)
    _M = Symbol(string("_", sym))
    @eval begin
        $M(iter::SubObjectIterator) = $_M(Val(iter.Acc), iter.Obj; iter.options...)
        $_M(::Any, ::Polymake.BigObject) = throw(ArgumentError(string($name, " not defined in this context.")))
    end
end

# Matrices with rational only elements
for (sym, name) in (("linear_inequality_matrix", "Linear Inequality Matrix"), ("affine_inequality_matrix", "Affine Inequality Matrix"), ("linear_equation_matrix", "Linear Equation Matrix"), ("affine_equation_matrix", "Affine Equation Matrix"))
    M = Symbol(sym)
    _M = Symbol(string("_", sym))
    @eval begin
        $M(iter::SubObjectIterator) = matrix(QQ, Matrix{fmpq}($_M(Val(iter.Acc), iter.Obj; iter.options...)))
        $_M(::Any, ::Polymake.BigObject) = throw(ArgumentError(string($name, " not defined in this context.")))
    end
end

# Matrices with rational or integer elements
for (sym, name) in (("point_matrix", "Point Matrix"), ("vector_matrix", "Vector Matrix"), ("generator_matrix", "Generator Matrix"))
    M = Symbol(sym)
    _M = Symbol(string("_", sym))
    @eval begin
        $M(iter::SubObjectIterator{<:AbstractVector{Polymake.Rational}}) = matrix(QQ, Matrix{fmpq}($_M(Val(iter.Acc), iter.Obj; iter.options...)))
        $M(iter::SubObjectIterator{<:AbstractVector{Polymake.Integer}}) = matrix(ZZ, $_M(Val(iter.Acc), iter.Obj; iter.options...))
        $_M(::Any, ::Polymake.BigObject) = throw(ArgumentError(string($name, " not defined in this context.")))
    end
end

function matrix_for_polymake(iter::SubObjectIterator; homogenized=false)
    if hasmethod(_matrix_for_polymake, Tuple{Val{iter.Acc}})
        return _matrix_for_polymake(Val(iter.Acc))(Val(iter.Acc), iter.Obj; homogenized=homogenized, iter.options...)
    else
        throw(ArgumentError("Matrix for Polymake not defined in this context."))
    end
end

# primitive generators only for ray based iterators
matrix(R::FlintIntegerRing, iter::SubObjectIterator{RayVector{Polymake.Rational}}) =
    matrix(R, Polymake.common.primitive(matrix_for_polymake(iter)))
matrix(R::FlintIntegerRing, iter::SubObjectIterator{<:Union{RayVector{Polymake.Integer},PointVector{Polymake.Integer}}}) =
    matrix(R, matrix_for_polymake(iter))
matrix(R::FlintRationalField, iter::SubObjectIterator{<:Union{RayVector,PointVector}}) =
    matrix(R, Matrix{fmpq}(matrix_for_polymake(iter)))

function linear_matrix_for_polymake(iter::SubObjectIterator)
    if hasmethod(_linear_matrix_for_polymake, Tuple{Val{iter.Acc}})
        return _linear_matrix_for_polymake(Val(iter.Acc))(Val(iter.Acc), iter.Obj; iter.options...)
    elseif hasmethod(_affine_matrix_for_polymake, Tuple{Val{iter.Acc}})
        res = _affine_matrix_for_polymake(Val(iter.Acc))(Val(iter.Acc), iter.Obj; iter.options...)
        iszero(res[:, 1]) || throw(ArgumentError("Input not linear."))
        return res[:, 2:end]
    end
    throw(ArgumentError("Linear Matrix for Polymake not defined in this context."))
end

function affine_matrix_for_polymake(iter::SubObjectIterator)
    if hasmethod(_affine_matrix_for_polymake, Tuple{Val{iter.Acc}})
        return _affine_matrix_for_polymake(Val(iter.Acc))(Val(iter.Acc), iter.Obj; iter.options...)
    elseif hasmethod(_linear_matrix_for_polymake, Tuple{Val{iter.Acc}})
        return homogenize(_linear_matrix_for_polymake(Val(iter.Acc))(Val(iter.Acc), iter.Obj; iter.options...), 0)
    end
    throw(ArgumentError("Affine Matrix for Polymake not defined in this context."))
end

function halfspace_matrix_pair(iter::SubObjectIterator)
    try
        h = affine_matrix_for_polymake(iter)
        return (A = matrix(QQ, Matrix{fmpq}(h[:, 2:end])), b = -h[:, 1])
    catch e
        throw(ArgumentError("Halfspace-Matrix-Pair not defined in this context."))
    end
end

Polymake.convert_to_pm_type(::Type{SubObjectIterator{RayVector{T}}}) where T = Polymake.Matrix{T}
Polymake.convert_to_pm_type(::Type{SubObjectIterator{PointVector{T}}}) where T = Polymake.Matrix{T}
Base.convert(::Type{<:Polymake.Matrix}, iter::SubObjectIterator) = matrix_for_polymake(iter; homogenized=true)

function homogenized_matrix(x::SubObjectIterator{<:PointVector}, v::Number = 1)
    if v != 1
        throw(ArgumentError("PointVectors can only be (re-)homogenized with parameter 1, please convert to a matrix first."))
    end
    return matrix_for_polymake(x; homogenized=true)
end
function homogenized_matrix(x::SubObjectIterator{<:RayVector}, v::Number = 0)
    if v != 0
        throw(ArgumentError("RayVectors can only be (re-)homogenized with parameter 0, please convert to a matrix first."))
    end
    return matrix_for_polymake(x; homogenized=true)
end
