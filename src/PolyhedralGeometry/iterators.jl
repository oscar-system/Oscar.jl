################################################################################
######## Scalar types
################################################################################
const scalar_types = Union{fmpq, nf_elem, Float64}

const scalar_type_to_oscar = Dict{String, Type}([("Rational", fmpq),
                                ("QuadraticExtension<Rational>", nf_elem),
                                ("Float", Float64)])

const scalar_type_to_polymake = Dict{Type, Type}([(fmpq, Polymake.Rational),
                                    (nf_elem, Polymake.QuadraticExtension{Polymake.Rational}),
                                    (Union{fmpq, nf_elem}, Polymake.QuadraticExtension{Polymake.Rational}),    # needed for Halfspace{nf_elem} etc
                                    (Float64, Float64)])

const scalar_types_extended = Union{scalar_types, fmpz}

################################################################################
######## Vector types
################################################################################

struct PointVector{U} <: AbstractVector{U}
    p::Vector{U}
    
    PointVector{U}(p::AbstractVector) where U<:scalar_types_extended = new{U}(p)
    PointVector(p::AbstractVector) = new{fmpq}(p)
end

Base.IndexStyle(::Type{<:PointVector}) = IndexLinear()

Base.getindex(po::PointVector{T}, i::Base.Integer) where T<:scalar_types_extended  = convert(T, po.p[i])

function Base.setindex!(po::PointVector, val, i::Base.Integer)
    @boundscheck checkbounds(po.p, i)
    po.p[i] = val
    return val
end

Base.firstindex(::PointVector) = 1
Base.lastindex(iter::PointVector) = length(iter)
Base.size(po::PointVector) = size(po.p)

PointVector{nf_elem}(p::AbstractVector) = PointVector{nf_scalar}(p)

PointVector{U}(n::Base.Integer) where U<:scalar_types_extended = PointVector{U}(zeros(U, n))

function Base.similar(X::PointVector, ::Type{S}, dims::Dims{1}) where S <: scalar_types_extended
    return PointVector{S}(dims...)
end

Base.BroadcastStyle(::Type{<:PointVector}) = Broadcast.ArrayStyle{PointVector}()

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{PointVector}}, ::Type{ElType}) where ElType
    return PointVector{ElType}(axes(bc)...)
end

################################################################################

struct RayVector{U} <: AbstractVector{U}
    p::Vector{U}
    
    RayVector{U}(p::AbstractVector) where U<:scalar_types_extended = new{U}(p)
    RayVector(p::AbstractVector) = new{fmpq}(p)
end

Base.IndexStyle(::Type{<:RayVector}) = IndexLinear()

Base.getindex(po::RayVector{T}, i::Base.Integer) where T<:scalar_types_extended  = convert(T, po.p[i])

function Base.setindex!(po::RayVector, val, i::Base.Integer)
    @boundscheck checkbounds(po.p, i)
    po.p[i] = val
    return val
end

Base.firstindex(::RayVector) = 1
Base.lastindex(iter::RayVector) = length(iter)
Base.size(po::RayVector) = size(po.p)

RayVector{nf_elem}(p::AbstractVector) = RayVector{nf_scalar}(p)

RayVector{U}(n::Base.Integer) where U<:scalar_types_extended = RayVector{U}(zeros(U, n))

function Base.similar(X::RayVector, ::Type{S}, dims::Dims{1}) where S <: scalar_types_extended
    return RayVector{S}(dims...)
end

Base.BroadcastStyle(::Type{<:RayVector}) = Broadcast.ArrayStyle{RayVector}()

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{RayVector}}, ::Type{ElType}) where ElType
    return RayVector{ElType}(axes(bc)...)
end

################################################################################

Base.:*(k::scalar_types_extended, po::Union{PointVector, RayVector}) = k .* po

################################################################################
######## Halfspaces and Hyperplanes
################################################################################

abstract type Halfspace{T} end

################################################################################

@doc Markdown.doc"""
    Halfspace(a, b)

One halfspace `H(a,b)` is given by a vector `a` and a value `b` such that
$$H(a,b) = \{ x | ax ≤ b \}.$$
"""
struct AffineHalfspace{T} <: Halfspace{T}
    a::Vector{T}
    b::T
    
    AffineHalfspace{T}(a::Union{MatElem, AbstractMatrix, AbstractVector}, b=0) where T<:scalar_types = new{T}(vec(a), b)
    AffineHalfspace(a::Union{MatElem, AbstractMatrix, AbstractVector}, b=0) = new{fmpq}(vec(a), b)
end

Halfspace(a, b) = AffineHalfspace(a, b)
Halfspace{T}(a, b) where T<:scalar_types = AffineHalfspace{T}(a, b)

invert(H::AffineHalfspace{T}) where T<:scalar_types = AffineHalfspace{T}(-normal_vector(H), -negbias(H))

################################################################################

struct LinearHalfspace{T} <: Halfspace{T}
    a::Vector{T}
    
    LinearHalfspace{T}(a::Union{MatElem, AbstractMatrix, AbstractVector}) where T<:scalar_types = new{T}(vec(a))
    LinearHalfspace(a::Union{MatElem, AbstractMatrix, AbstractVector}) = new{fmpq}(vec(a))
end

Halfspace(a) = LinearHalfspace(a)
Halfspace{T}(a) where T<:scalar_types = LinearHalfspace{T}(a)

invert(H::LinearHalfspace{T}) where T<:scalar_types = LinearHalfspace{T}(-normal_vector(H))

################################################################################

abstract type Hyperplane{T} end

################################################################################

@doc Markdown.doc"""
    AffineHyperplane(a, b)

One hyperplane `H(a,b)` is given by a vector `a` and a value `b` such that
$$H(a,b) = \{ x | ax = b \}.$$
"""
struct AffineHyperplane{T} <: Hyperplane{T}
    a::Vector{T}
    b::T
    
    AffineHyperplane{T}(a::Union{MatElem, AbstractMatrix, AbstractVector}, b=0) where T<:scalar_types = new{T}(vec(a), b)
    AffineHyperplane(a::Union{MatElem, AbstractMatrix, AbstractVector}, b=0) = new{fmpq}(vec(a), b)
end

Hyperplane(a, b) = AffineHyperplane(a, b)
Hyperplane{T}(a, b) where T<:scalar_types = AffineHyperplane{T}(a, b)

################################################################################

struct LinearHyperplane{T} <: Hyperplane{T}
    a::Vector{T}
    
    LinearHyperplane{T}(a::Union{MatElem, AbstractMatrix, AbstractVector}) where T<:scalar_types = new{T}(vec(a))
    LinearHyperplane(a::Union{MatElem, AbstractMatrix, AbstractVector}) = new{fmpq}(vec(a))
end

Hyperplane(a) = LinearHyperplane(a)
Hyperplane{T}(a) where T<:scalar_types = LinearHyperplane{T}(a)

################################################################################

#  Field access
negbias(H::Union{AffineHalfspace{T}, AffineHyperplane{T}}) where T<:scalar_types = H.b
negbias(H::Union{LinearHalfspace{T}, LinearHyperplane{T}}) where T<:scalar_types = T(0)
normal_vector(H::Union{Halfspace{T}, Hyperplane{T}}) where T <: scalar_types = Vector{T}(H.a)

ambient_dim(x::Union{Halfspace, Hyperplane}) = length(x.a)

# TODO: abstract notion of equality
Base.:(==)(x::AffineHalfspace, y::AffineHalfspace) = x.a == y.a && x.b == y.b

Base.:(==)(x::LinearHalfspace, y::LinearHalfspace) = x.a == y.a

Base.:(==)(x::AffineHyperplane, y::AffineHyperplane) = x.a == y.a && x.b == y.b

Base.:(==)(x::LinearHyperplane, y::LinearHyperplane) = x.a == y.a

################################################################################

for T in [LinearHalfspace, LinearHyperplane]
    @eval begin
        $T{nf_elem}(a::Union{MatElem, AbstractMatrix, AbstractVector}) = $T{nf_scalar}(a)
    end
end

for T in [AffineHalfspace, AffineHyperplane]
    @eval begin
        $T{nf_elem}(a::Union{MatElem, AbstractMatrix, AbstractVector}, b = 0) = $T{nf_scalar}(a, b)
    end
end

################################################################################
######## SubObjectIterator
################################################################################

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

# `options` is empty by default
SubObjectIterator{T}(Obj::Polymake.BigObject, Acc::Function, n::Base.Integer) where T = SubObjectIterator{T}(Obj, Acc, n, NamedTuple())

# Force `nf_scalar` for this `Pair` descpription of `Halfspace`s/`Hyperplane`s
# derived from `nf_elem`templated object
SubObjectIterator{Pair{Matrix{nf_elem}, nf_elem}}(Obj::Polymake.BigObject, Acc::Function, n::Base.Integer, options::NamedTuple = NamedTuple()) = SubObjectIterator{Pair{Matrix{nf_scalar}, nf_scalar}}(Obj, Acc, n, options)

Base.IndexStyle(::Type{<:SubObjectIterator}) = IndexLinear()

function Base.getindex(iter::SubObjectIterator{T}, i::Base.Integer) where T
    @boundscheck 1 <= i && i <= iter.n
    return iter.Acc(T, iter.Obj, i; iter.options...)
end

Base.firstindex(::SubObjectIterator) = 1
Base.lastindex(iter::SubObjectIterator) = length(iter)
Base.size(iter::SubObjectIterator) = (iter.n,)

################################################################################

# Incidence matrices
for (sym, name) in (("ray_indices", "Incidence Matrix resp. rays"), ("vertex_indices", "Incidence Matrix resp. vertices"), ("vertex_and_ray_indices", "Incidence Matrix resp. vertices and rays"))
    M = Symbol(sym)
    _M = Symbol(string("_", sym))
    @eval begin
        $M(iter::SubObjectIterator) = $_M(Val(iter.Acc), iter.Obj; iter.options...)
        $_M(::Any, ::Polymake.BigObject) = throw(ArgumentError(string($name, " not defined in this context.")))
    end
end

# Matrices with rational or integer elements
for (sym, name) in (("point_matrix", "Point Matrix"), ("vector_matrix", "Vector Matrix"), ("generator_matrix", "Generator Matrix"))
    M = Symbol(sym)
    _M = Symbol(string("_", sym))
    @eval begin
        $M(iter::SubObjectIterator{<:AbstractVector{fmpq}}) = matrix(QQ, Matrix{fmpq}($_M(Val(iter.Acc), iter.Obj; iter.options...)))
        $M(iter::SubObjectIterator{<:AbstractVector{fmpz}}) = matrix(ZZ, $_M(Val(iter.Acc), iter.Obj; iter.options...))
        $M(iter::SubObjectIterator{<:AbstractVector{nf_elem}}) = Matrix{nf_scalar}($_M(Val(iter.Acc), iter.Obj; iter.options...))
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
matrix(R::FlintIntegerRing, iter::SubObjectIterator{RayVector{fmpq}}) =
    matrix(R, Polymake.common.primitive(matrix_for_polymake(iter)))
matrix(R::FlintIntegerRing, iter::SubObjectIterator{<:Union{RayVector{fmpz},PointVector{fmpz}}}) =
    matrix(R, matrix_for_polymake(iter))
matrix(R::FlintRationalField, iter::SubObjectIterator{<:Union{RayVector{fmpq}, PointVector{fmpq}}}) =
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

homogenized_matrix(::SubObjectIterator, v::Number) = throw(ArgumentError("Content of SubObjectIterator not suitable for homogenized_matrix."))

unhomogenized_matrix(x::SubObjectIterator{<:RayVector}) = matrix_for_polymake(x)

unhomogenized_matrix(x::SubObjectIterator{<:PointVector}) = throw(ArgumentError("unhomogenized_matrix only meaningful for RayVectors"))

ambient_dim(x::SubObjectIterator) = Polymake.polytope.ambient_dim(x.Obj)

################################################################################
######## Unify matrices
################################################################################

# homogenized_matrix -> matrix_for_polymake
# linear_matrix_for_polymake
const SomeMatrix = Union{AnyVecOrMat, SubObjectIterator}

# affine_matrix_for_polymake
const LinearExpressionSet = Union{SubObjectIterator, Tuple{SomeMatrix, Any}}
