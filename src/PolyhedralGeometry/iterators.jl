################################################################################
######## Vector types
################################################################################

struct PointVector{U} <: AbstractVector{U}
    p::Vector{U}
    
    PointVector{U}(p::AbstractVector) where U<:scalar_types_extended = new{U}(p)
    PointVector(p::AbstractVector) = new{QQFieldElem}(p)
end

Base.IndexStyle(::Type{<:PointVector}) = IndexLinear()

Base.getindex(po::PointVector{T}, i::Base.Integer) where T<:scalar_types_extended = po.p[i]

function Base.setindex!(po::PointVector, val, i::Base.Integer)
    @boundscheck checkbounds(po.p, i)
    po.p[i] = val
    return val
end

Base.firstindex(::PointVector) = 1
Base.lastindex(iter::PointVector) = length(iter)
Base.size(po::PointVector) = size(po.p)

# Forward multiplication with oscar matrices.
Base.:*(A::MatElem, v::PointVector) = A*v.p

PointVector{U}(p::Union{Field, ZZRing}, v::AbstractVector) where U<:scalar_types_extended = PointVector{U}(p.(v))

PointVector{U}(p::Union{Field, ZZRing}, n::Base.Integer) where U<:scalar_types_extended = PointVector{U}(p.(zeros(Int, n)))

PointVector{U}(::UndefInitializer, n::Base.Integer) where U<:scalar_types_extended = PointVector{U}(Vector{U}(undef, n))

PointVector{U}(n::UR) where {U<:scalar_types_extended, UR<:Base.AbstractUnitRange} = PointVector{U}(undef, length(n))

function Base.similar(X::PointVector, ::Type{S}, dims::Dims{1}) where S <: scalar_types_extended
    return PointVector{S}(undef, dims...)
end

Base.BroadcastStyle(::Type{<:PointVector}) = Broadcast.ArrayStyle{PointVector}()

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{PointVector}}, ::Type{ElType}) where ElType
    return PointVector{ElType}(axes(bc)...)
end

################################################################################

struct RayVector{U} <: AbstractVector{U}
    p::Vector{U}
    
    RayVector{U}(p::AbstractVector) where U<:scalar_types_extended = new{U}(p)
    RayVector(p::AbstractVector) = new{QQFieldElem}(p)
end

Base.IndexStyle(::Type{<:RayVector}) = IndexLinear()

Base.getindex(po::RayVector{T}, i::Base.Integer) where T<:scalar_types_extended = po.p[i]

function Base.setindex!(po::RayVector, val, i::Base.Integer)
    @boundscheck checkbounds(po.p, i)
    po.p[i] = val
    return val
end

Base.firstindex(::RayVector) = 1
Base.lastindex(iter::RayVector) = length(iter)
Base.size(po::RayVector) = size(po.p)

# Forward multiplication with oscar matrices.
Base.:*(A::MatElem, v::RayVector) = A*v.p

RayVector{U}(p::Union{Field, ZZRing}, v::AbstractVector) where U<:scalar_types_extended = RayVector{U}(p.(v))

RayVector{U}(p::Union{Field, ZZRing}, n::Base.Integer) where U<:scalar_types_extended = RayVector{U}(p.(zeros(Int, n)))

RayVector{U}(::UndefInitializer, n::Base.Integer) where U<:scalar_types_extended = RayVector{U}(Vector{U}(undef, n))

RayVector{U}(n::UR) where {U<:scalar_types_extended, UR<:Base.AbstractUnitRange} = RayVector{U}(undef, length(n))

function Base.similar(X::RayVector, ::Type{S}, dims::Dims{1}) where S <: scalar_types_extended
    return RayVector{S}(undef, dims...)
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

@doc raw"""
    Halfspace(a, b)

One halfspace `H(a,b)` is given by a vector `a` and a value `b` such that
$$H(a,b) = \{ x | ax ≤ b \}.$$
"""
struct AffineHalfspace{T} <: Halfspace{T}
    a::Vector{T}
    b::T
    parent_field::Field

    AffineHalfspace{T}(p::Field, a::Union{MatElem, AbstractMatrix, AbstractVector}, b) where T<:scalar_types = new{T}(p.(vec(a)), p(b), p)
end

halfspace(a::Union{MatElem, AbstractMatrix, AbstractVector}, b) = affine_halfspace(a, b)
halfspace(f::Union{Type{T}, Field}, a::Union{MatElem, AbstractMatrix, AbstractVector}, b) where T<:scalar_types = affine_halfspace(f, a, b)

invert(H::AffineHalfspace{T}) where T<:scalar_types = AffineHalfspace{T}(coefficient_field(H), -normal_vector(H), -negbias(H))

function affine_halfspace(f::Union{Type{T}, Field}, a::Union{MatElem, AbstractMatrix, AbstractVector}, b = 0) where T<:scalar_types
    parent_field, scalar_type = _determine_parent_and_scalar(f, a, b)
    return AffineHalfspace{scalar_type}(parent_field, a, b)
end

affine_halfspace(a::Union{MatElem, AbstractMatrix, AbstractVector}, b = 0) = affine_halfspace(QQ, a, b)

################################################################################

struct LinearHalfspace{T} <: Halfspace{T}
    a::Vector{T}
    parent_field::Field
    
    LinearHalfspace{T}(p::Field, a::Union{MatElem, AbstractMatrix, AbstractVector}) where T<:scalar_types = new{T}(p.(vec(a)), p)
end

halfspace(a::Union{MatElem, AbstractMatrix, AbstractVector}) = linear_halfspace(a)
halfspace(f::Union{Type{T}, Field}, a::Union{MatElem, AbstractMatrix, AbstractVector}) where T<:scalar_types = linear_halfspace(f, a)

invert(H::LinearHalfspace{T}) where T<:scalar_types = LinearHalfspace{T}(coefficient_field(H), -normal_vector(H))

function linear_halfspace(f::Union{Type{T}, Field}, a::Union{MatElem, AbstractMatrix, AbstractVector}) where T<:scalar_types
    parent_field, scalar_type = _determine_parent_and_scalar(f, a)
    return LinearHalfspace{scalar_type}(parent_field, a)
end

linear_halfspace(a::Union{MatElem, AbstractMatrix, AbstractVector}) = linear_halfspace(QQ, a)

coefficient_field(h::Halfspace) = h.parent_field

################################################################################

abstract type Hyperplane{T} end

################################################################################

@doc raw"""
    AffineHyperplane(a, b)

One hyperplane `H(a,b)` is given by a vector `a` and a value `b` such that
$$H(a,b) = \{ x | ax = b \}.$$
"""
struct AffineHyperplane{T} <: Hyperplane{T}
    a::Vector{T}
    b::T
    parent_field::Field
    
    AffineHyperplane{T}(p::Field, a::Union{MatElem, AbstractMatrix, AbstractVector}, b=0) where T<:scalar_types = new{T}(p.(vec(a)), p(b), p)
end

hyperplane(a::Union{MatElem, AbstractMatrix, AbstractVector}, b) = affine_hyperplane(a, b)
hyperplane(f::Union{Type{T}, Field}, a::Union{MatElem, AbstractMatrix, AbstractVector}, b) where T<:scalar_types = affine_hyperplane(f, a, b)

function affine_hyperplane(f::Union{Type{T}, Field}, a::Union{MatElem, AbstractMatrix, AbstractVector}, b = 0) where T<:scalar_types
    parent_field, scalar_type = _determine_parent_and_scalar(f, a, b)
    return AffineHyperplane{scalar_type}(parent_field, a, b)
end

affine_hyperplane(a::Union{MatElem, AbstractMatrix, AbstractVector}, b = 0) = affine_hyperplane(QQ, a, b)

################################################################################

struct LinearHyperplane{T} <: Hyperplane{T}
    a::Vector{T}
    parent_field::Field
    
    LinearHyperplane{T}(p::Field, a::Union{MatElem, AbstractMatrix, AbstractVector}) where T<:scalar_types = new{T}(p.(vec(a)), p)
end

hyperplane(a::Union{MatElem, AbstractMatrix, AbstractVector}) = linear_hyperplane(a)
hyperplane(f::Union{Type{T}, Field}, a::Union{MatElem, AbstractMatrix, AbstractVector}) where T<:scalar_types = linear_hyperplane(f, a)

function linear_hyperplane(f::Union{Type{T}, Field}, a::Union{MatElem, AbstractMatrix, AbstractVector}) where T<:scalar_types
    parent_field, scalar_type = _determine_parent_and_scalar(f, a)
    return LinearHyperplane{scalar_type}(parent_field, a)
end

linear_hyperplane(a::Union{MatElem, AbstractMatrix, AbstractVector}) = linear_hyperplane(QQ, a)

coefficient_field(h::Hyperplane) = h.parent_field

################################################################################

#  Field access
negbias(H::Union{AffineHalfspace{T}, AffineHyperplane{T}}) where T<:scalar_types = H.b
negbias(H::Union{LinearHalfspace{T}, LinearHyperplane{T}}) where T<:scalar_types = coefficient_field(H)(0)
normal_vector(H::Union{Halfspace{T}, Hyperplane{T}}) where T <: scalar_types = Vector{T}(H.a)

_ambient_dim(x::Union{Halfspace, Hyperplane}) = length(x.a)

# TODO: abstract notion of equality
Base.:(==)(x::AffineHalfspace, y::AffineHalfspace) = x.a == y.a && x.b == y.b

Base.:(==)(x::LinearHalfspace, y::LinearHalfspace) = x.a == y.a

Base.:(==)(x::AffineHyperplane, y::AffineHyperplane) = x.a == y.a && x.b == y.b

Base.:(==)(x::LinearHyperplane, y::LinearHyperplane) = x.a == y.a

Base.hash(x::T, h::UInt) where {T<:Union{AffineHalfspace,AffineHyperplane}} =
  hash((x.a, x.b), hash(T, h))

Base.hash(x::T, h::UInt) where {T<:Union{LinearHalfspace,LinearHyperplane}} =
  hash(x.a, hash(T, h))

################################################################################
######## SubObjectIterator
################################################################################

@doc raw"""
    SubObjectIterator(Obj, Acc, n, [options])

An iterator over a designated property of an object `Obj::PolyhedralObject` from Polyhedral Geometry.

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
    Obj::PolyhedralObjectUnion
    Acc::Function
    n::Int
    options::NamedTuple
end

# `options` is empty by default
SubObjectIterator{T}(Obj::PolyhedralObjectUnion, Acc::Function, n::Base.Integer) where T = SubObjectIterator{T}(Obj, Acc, n, NamedTuple())

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
for (sym, name) in (("facet_indices", "Incidence matrix resp. facets"), ("ray_indices", "Incidence Matrix resp. rays"), ("vertex_indices", "Incidence Matrix resp. vertices"), ("vertex_and_ray_indices", "Incidence Matrix resp. vertices and rays"))
    M = Symbol(sym)
    _M = Symbol("_", sym)
    @eval begin
        $M(iter::SubObjectIterator) = $_M(Val(iter.Acc), iter.Obj; iter.options...)
        $_M(::Any, ::PolyhedralObjectUnion) = throw(ArgumentError(string($name, " not defined in this context.")))
    end
end

# Matrices with rational or integer elements
for (sym, name) in (("point_matrix", "Point Matrix"), ("vector_matrix", "Vector Matrix"), ("generator_matrix", "Generator Matrix"))
    M = Symbol(sym)
    _M = Symbol("_", sym)
    @eval begin
        $M(iter::SubObjectIterator{<:AbstractVector{QQFieldElem}}) = matrix(QQ, $_M(Val(iter.Acc), iter.Obj; iter.options...))
        $M(iter::SubObjectIterator{<:AbstractVector{ZZRingElem}}) = matrix(ZZ, $_M(Val(iter.Acc), iter.Obj; iter.options...))
        $M(iter::SubObjectIterator{<:AbstractVector{<:FieldElem}}) = matrix(coefficient_field(iter.Obj), $_M(Val(iter.Acc), iter.Obj; iter.options...))
        $_M(::Any, ::PolyhedralObjectUnion) = throw(ArgumentError(string($name, " not defined in this context.")))
    end
end

function matrix_for_polymake(iter::SubObjectIterator; homogenized=false)
    if hasmethod(_matrix_for_polymake, Tuple{Val{iter.Acc}})
        return _matrix_for_polymake(Val(iter.Acc))(Val(iter.Acc), iter.Obj; homogenized=homogenized, iter.options...)
    else
        throw(ArgumentError("Matrix for Polymake not defined in this context."))
    end
end

function IncidenceMatrix(iter::SubObjectIterator)
    if hasmethod(_incidencematrix, Tuple{Val{iter.Acc}})
        return _incidencematrix(Val(iter.Acc))(Val(iter.Acc), iter.Obj; iter.options...)
    else
        throw(ArgumentError("IncidenceMatrix not defined in this context."))
    end
end

# primitive generators only for ray based iterators
matrix(R::ZZRing, iter::SubObjectIterator{RayVector{QQFieldElem}}) =
    matrix(R, Polymake.common.primitive(matrix_for_polymake(iter)))
matrix(R::ZZRing, iter::SubObjectIterator{<:Union{RayVector{ZZRingElem},PointVector{ZZRingElem}}}) =
    matrix(R, matrix_for_polymake(iter))
matrix(R::QQField, iter::SubObjectIterator{<:Union{RayVector{QQFieldElem}, PointVector{QQFieldElem}}}) =
    matrix(R, matrix_for_polymake(iter))
matrix(K, iter::SubObjectIterator{<:Union{RayVector{<:FieldElem}, PointVector{<:FieldElem}}}) =
    matrix(K, matrix_for_polymake(iter))


function linear_matrix_for_polymake(iter::SubObjectIterator)
    if hasmethod(_linear_matrix_for_polymake, Tuple{Val{iter.Acc}})
        return _linear_matrix_for_polymake(Val(iter.Acc))(Val(iter.Acc), iter.Obj; iter.options...)
    elseif hasmethod(_affine_matrix_for_polymake, Tuple{Val{iter.Acc}})
        res = _affine_matrix_for_polymake(Val(iter.Acc))(Val(iter.Acc), iter.Obj; iter.options...)
        @req iszero(res[:, 1]) "Input not linear."
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
Base.convert(::Type{<:Polymake.Matrix}, iter::SubObjectIterator) = assure_matrix_polymake(matrix_for_polymake(iter; homogenized=true))

function homogenized_matrix(x::SubObjectIterator{<:PointVector}, v::Number = 1)
    @req v == 1 "PointVectors can only be (re-)homogenized with parameter 1, please convert to a matrix first"
    return matrix_for_polymake(x; homogenized=true)
end
function homogenized_matrix(x::SubObjectIterator{<:RayVector}, v::Number = 0)
    @req v == 0 "RayVectors can only be (re-)homogenized with parameter 0, please convert to a matrix first"
    return matrix_for_polymake(x; homogenized=true)
end

function homogenized_matrix(x::AbstractVector{<:PointVector}, v::Number = 1)
    @req v == 1 "PointVectors can only be (re-)homogenized with parameter 1, please convert to a matrix first"
    return stack((homogenize(x[i], v) for i in 1:length(x))...)
end
function homogenized_matrix(x::AbstractVector{<:RayVector}, v::Number = 0)
    @req v == 0 "RayVectors can only be (re-)homogenized with parameter 0, please convert to a matrix first"
    return stack((homogenize(x[i], v) for i in 1:length(x))...)
end

homogenized_matrix(::SubObjectIterator, v::Number) = throw(ArgumentError("Content of SubObjectIterator not suitable for homogenized_matrix."))

unhomogenized_matrix(x::SubObjectIterator{<:RayVector}) = matrix_for_polymake(x)

unhomogenized_matrix(x::AbstractVector{<:PointVector}) = throw(ArgumentError("unhomogenized_matrix only meaningful for RayVectors"))

_ambient_dim(x::SubObjectIterator) = Polymake.polytope.ambient_dim(pm_object(x.Obj))

################################################################################

# Lineality often causes certain collections to be empty;
# the following definition allows to easily construct a working empty SOI

_empty_access() = nothing

function _empty_subobjectiterator(::Type{T}, Obj::PolyhedralObjectUnion) where T
    return SubObjectIterator{T}(Obj, _empty_access, 0, NamedTuple())
end

for f in ("_point_matrix", "_vector_matrix", "_generator_matrix")
    M = Symbol(f)
    @eval begin
        function $M(::Val{_empty_access}, P::PolyhedralObjectUnion; homogenized=false)
            scalar_regexp = match(r"[^<]*<(.*)>[^>]*", String(Polymake.type_name(pm_object(P))))
            typename = scalar_regexp[1]
            T = _scalar_type_to_polymake(scalar_type_to_oscar[typename])
            return Polymake.Matrix{T}(undef, 0, Polymake.polytope.ambient_dim(pm_object(P)) + homogenized)
        end
    end
end

for f in ("_facet_indices", "_ray_indices", "_vertex_indices", "_vertex_and_ray_indices")
    M = Symbol(f)
    @eval begin
        $M(::Val{_empty_access}, P::PolyhedralObjectUnion) = return Polymake.IncidenceMatrix(0, Polymake.polytope.ambient_dim(P))
    end
end

for f in ("_linear_inequality_matrix", "_linear_equation_matrix")
    M = Symbol(f)
    @eval begin
        function $M(::Val{_empty_access}, P::PolyhedralObjectUnion)
            scalar_regexp = match(r"[^<]*<(.*)>[^>]*", String(Polymake.type_name(P)))
            typename = scalar_regexp[1]
            T = _scalar_type_to_polymake(scalar_type_to_oscar[typename])
            return Polymake.Matrix{T}(undef, 0, Polymake.polytope.ambient_dim(P))
        end
    end
end

for f in ("_affine_inequality_matrix", "_affine_equation_matrix")
    M = Symbol(f)
    @eval begin
        function $M(::Val{_empty_access}, P::PolyhedralObjectUnion)
            scalar_regexp = match(r"[^<]*<(.*)>[^>]*", String(Polymake.type_name(P)))
            typename = scalar_regexp[1]
            T = _scalar_type_to_polymake(scalar_type_to_oscar[typename])
            return Polymake.Matrix{T}(undef, 0, Polymake.polytope.ambient_dim(P) + 1)
        end
    end
end

_matrix_for_polymake(::Val{_empty_access}) = _point_matrix

################################################################################
######## Unify matrices
################################################################################

# vector-like types: (un-)homogenized_matrix -> matrix_for_polymake
# linear types: linear_matrix_for_polymake 
# affine types: affine_matrix_for_polymake
const AbstractCollection = Dict{UnionAll, Union}([(PointVector, AnyVecOrMat),
                                                    (RayVector, AnyVecOrMat),
                                                    (LinearHalfspace, Union{AbstractVector{<:Halfspace}, SubObjectIterator{<:Halfspace}, AnyVecOrMat}),
                                                    (LinearHyperplane, Union{AbstractVector{<:Hyperplane}, SubObjectIterator{<:Hyperplane}, AnyVecOrMat}),
                                                    (AffineHalfspace, Union{AbstractVector{<:Halfspace}, SubObjectIterator{<:Halfspace}, Tuple{<:AnyVecOrMat, <:Any}}),
                                                    (AffineHyperplane, Union{AbstractVector{<:Hyperplane}, SubObjectIterator{<:Hyperplane}, Tuple{<:AnyVecOrMat, <:Any}})])

affine_matrix_for_polymake(x::Union{Halfspace, Hyperplane}) = stack(augment(normal_vector(x), -negbias(x)))
affine_matrix_for_polymake(x::AbstractVector{<:Union{Halfspace, Hyperplane}}) = stack((affine_matrix_for_polymake(x[i]) for i in 1:length(x))...)
linear_matrix_for_polymake(x::Union{Halfspace, Hyperplane}) = negbias(x) == 0 ? stack(normal_vector(x)) : throw(ArgumentError("Input not linear."))
linear_matrix_for_polymake(x::AbstractVector{<:Union{Halfspace, Hyperplane}}) = stack((linear_matrix_for_polymake(x[i]) for i in 1:length(x))...)
