################################################################################
######## Vector types
################################################################################

for (T, _t) in ((:PointVector, :point_vector), (:RayVector, :ray_vector))
  @eval begin
    struct $T{U} <: AbstractVector{U}
      p::MatElem{U}

      $T{U}(p::MatElem{U}) where {U<:scalar_types_extended} = new{U}(p)
    end

    Base.IndexStyle(::Type{<:$T}) = IndexLinear()
    Base.getindex(po::$T{U}, i::Base.Integer) where {U} = po.p[1, i]::U

    function Base.setindex!(po::$T, val, i::Base.Integer)
      @boundscheck 1 <= length(po) <= i
      po.p[1, i] = val
      return val
    end

    Base.firstindex(::$T) = 1
    Base.lastindex(iter::$T) = length(iter)
    Base.size(po::$T) = (size(po.p, 2)::Int,)

    coefficient_field(po::$T) = base_ring(po.p)

    function $_t(p::Union{scalar_type_or_field,ZZRing}, v::AbstractVector)
      parent_field, scalar_type = _determine_parent_and_scalar(p, v)
      n = length(v)
      mat = matrix(parent_field, 1, n, collect(v)) # collect: workaround for constructor typing
      return $T{scalar_type}(mat)
    end

    function $_t(p::Union{scalar_type_or_field,ZZRing}, n::Base.Integer)
      parent_field, scalar_type = _determine_parent_and_scalar(p)
      mat = zero_matrix(parent_field, 1, n)
      return $T{scalar_type}(mat)
    end

    $_t(x::Union{AbstractVector,Base.Integer}) = $_t(QQ, x)

    function Base.similar(X::$T, ::Type{S}, dims::Dims{1}) where {S<:scalar_types_extended}
      return $_t(coefficient_field(X), dims...)
    end

    Base.BroadcastStyle(::Type{<:$T}) = Broadcast.ArrayStyle{$T}()

    _parent_or_coefficient_field(::Type{TT}, po::$T{<:TT}) where {TT<:FieldElem} =
      coefficient_field(po)
    _find_elem_type(po::$T) = elem_type(coefficient_field(po))

    function Base.similar(
      bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{$T}},
      ::Type{<:Union{scalar_types,ZZRingElem}})
      e = bc.f(first.(bc.args)...)
      return $_t(parent(e), axes(bc)...)
    end

    function Base.similar(
      bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{$T}}, ::Type{ElType}
    ) where {ElType}
      return Vector{ElType}(undef, length(axes(bc)...))
    end

    Base.:*(k::scalar_types_extended, po::$T) = k .* po

    Base.:*(A::MatElem, v::$T) = A * transpose(v.p)
  end
end

@doc """
    point_vector(p = QQ, v::AbstractVector)

Return a `PointVector` resembling a point whose coordinates equal the entries of `v`.
`p` specifies the `Field` or `Type` of its coefficients.
"""
point_vector

@doc """
    ray_vector(p = QQ, v::AbstractVector)

Return a `RayVector` resembling a ray from the origin through the point whose coordinates equal the entries of `v`.
`p` specifies the `Field` or `Type` of its coefficients.
"""
ray_vector

function Base.:(==)(x::RayVector, y::RayVector)
  ix = findfirst(!is_zero, x)
  iy = findfirst(!is_zero, y)
  ix == iy || return false
  isnothing(ix) && return true
  sign(x[ix]) == sign(y[iy]) || return false
  return y[iy] * x.p == x[ix] * y.p
end

function Base.:(==)(x::RayVector, y::AbstractVector)
  ry = ray_vector(coefficient_field(x), y)
  return x == ry
end

Base.:(==)(x::AbstractVector, y::RayVector) = y == x

Base.:(==)(::PointVector, ::RayVector) =
  throw(ArgumentError("Cannot compare PointVector to RayVector"))
Base.:(==)(::RayVector, ::PointVector) =
  throw(ArgumentError("Cannot compare PointVector to RayVector"))

Base.isequal(x::RayVector, y::RayVector) = x == y

Base.hash(x::RayVector, h::UInt) = hash(collect(sign.(x)), hash(coefficient_field(x), h))

################################################################################
######## Halfspaces and Hyperplanes
################################################################################

for (h, comp) in (("halfspace", "≤"), ("hyperplane", "="))
  H = uppercasefirst(h)
  Habs = Symbol(H)
  Haff = Symbol("Affine", H)
  Hlin = Symbol("Linear", H)
  Fabs = Symbol(h)
  Faff = Symbol("affine_", h)
  Flin = Symbol("linear_", h)

  @eval begin
    abstract type $Habs{T} end

    # Affine types
    struct $Haff{T} <: $Habs{T}
      a::MatElem{T}
      b::T

      $Haff{T}(a::MatElem{T}, b::T) where {T<:scalar_types} = new{T}(a, b)
    end

    $Fabs(a::Union{MatElem,AbstractMatrix,AbstractVector}, b) = $Faff(a, b)
    $Fabs(f::scalar_type_or_field, a::Union{MatElem,AbstractMatrix,AbstractVector}, b) =
      $Faff(f, a, b)

    invert(H::$Haff{T}) where {T<:scalar_types} = $Haff{T}(-H.a, -negbias(H))

    @doc """
    $($Faff)(p = QQ, a, b)

Return the `$($Haff)` `H(a,b)`, which is given by a vector `a` and a value `b` such that
\$\$H(a,b) = \\{ x | ax $($comp) b \\}.\$\$
`p` specifies the `Field` or `Type` of its coefficients.
    """
    function $Faff(
      f::scalar_type_or_field, a::Union{MatElem,AbstractMatrix,AbstractVector}, b=0
    )
      parent_field, scalar_type = _determine_parent_and_scalar(f, a, b)
      mat = matrix(parent_field, 1, length(a), collect(a))
      return $Haff{scalar_type}(mat, parent_field(b))
    end

    $Faff(a::Union{MatElem,AbstractMatrix,AbstractVector}, b=0) = $Faff(QQ, a, b)

    # Linear types

    struct $Hlin{T} <: $Habs{T}
      a::MatElem{T}

      $Hlin{T}(a::MatElem{T}) where {T<:scalar_types} = new{T}(a)
    end

    $Fabs(a::Union{MatElem,AbstractMatrix,AbstractVector}) = $Flin(a)
    $Fabs(f::scalar_type_or_field, a::Union{MatElem,AbstractMatrix,AbstractVector}) =
      $Flin(f, a)

    invert(H::$Hlin{T}) where {T<:scalar_types} = $Hlin{T}(-H.a)

    @doc """
    $($Flin)(p = QQ, a, b)

Return the `$($Hlin)` `H(a)`, which is given by a vector `a` such that
\$\$H(a,b) = \\{ x | ax $($comp) 0 \\}.\$\$
`p` specifies the `Field` or `Type` of its coefficients.
    """
    function $Flin(f::scalar_type_or_field, a::Union{MatElem,AbstractMatrix,AbstractVector})
      parent_field, scalar_type = _determine_parent_and_scalar(f, a)
      mat = matrix(parent_field, 1, length(a), collect(a))
      return $Hlin{scalar_type}(mat)
    end

    $Flin(a::Union{MatElem,AbstractMatrix,AbstractVector}) = $Flin(QQ, a)

    coefficient_field(h::$Habs) = base_ring(h.a)

    _find_elem_type(h::$Habs) = elem_type(coefficient_field(h))
    _parent_or_coefficient_field(::Type{T}, h::$Habs{<:T}) where {T<:FieldElem} =
      coefficient_field(h)
  end
end

#  Field access
negbias(H::Union{AffineHalfspace,AffineHyperplane}) = H.b
negbias(H::Union{LinearHalfspace,LinearHyperplane}) = coefficient_field(H)(0)
normal_vector(H::Union{Halfspace,Hyperplane}) = H.a[1, :]

_ambient_dim(x::Union{Halfspace,Hyperplane}) = length(x.a)

function Base.:(==)(x::Halfspace, y::Halfspace)
  ax = normal_vector(x)
  ay = normal_vector(y)
  ix = findfirst(!is_zero, ax)
  iy = findfirst(!is_zero, ay)
  ix == iy || return false
  r = y.a[iy]//x.a[ix]
  r > 0 || return false
  return (r .* ax == ay) && (r * negbias(x) == negbias(y))
end

function Base.:(==)(x::Hyperplane, y::Hyperplane)
  ax = normal_vector(x)
  ay = normal_vector(y)
  ix = findfirst(!is_zero, ax)
  iy = findfirst(!is_zero, ay)
  ix == iy || return false
  r = y.a[iy]//x.a[ix]
  return (r .* ax == ay) && (r * negbias(x) == negbias(y))
end

Base.in(x::AbstractVector, y::Hyperplane) = (dot(x, normal_vector(y)) == negbias(y))
Base.in(x::AbstractVector, y::Halfspace) = (dot(x, normal_vector(y)) <= negbias(y))
# A ray vector needs a base point for containment in an affine space, so we
# just error when this combination is tested.
Base.in(x::RayVector, y::T) where {T<:Union{AffineHalfspace,
    AffineHyperplane}} =
  throw(ArgumentError("Containment of RayVector in affine spaces is not
                      well-defined."))

function Base.hash(x::Union{<:Halfspace,<:Hyperplane}, h::UInt)
  a = normal_vector(x)
  f = inv(first(Iterators.filter(!is_zero, a)))
  b = negbias(x)
  if x isa Halfspace
    f = abs(f)
  end
  hash((f .* a, f * b), h)
end

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
SubObjectIterator{T}(Obj::PolyhedralObjectUnion, Acc::Function, n::Base.Integer) where {T} =
  SubObjectIterator{T}(Obj, Acc, n, NamedTuple())

Base.IndexStyle(::Type{<:SubObjectIterator}) = IndexLinear()

function Base.getindex(iter::SubObjectIterator{T}, i::Base.Integer) where {T}
  @boundscheck 1 <= i && i <= iter.n
  return iter.Acc(T, iter.Obj, i; iter.options...)::T
end

Base.firstindex(::SubObjectIterator) = 1
Base.lastindex(iter::SubObjectIterator) = length(iter)
Base.size(iter::SubObjectIterator) = (iter.n,)

################################################################################

# Incidence matrices
for (sym, name) in (
  ("facet_indices", "Incidence matrix resp. facets"),
  ("ray_indices", "Incidence Matrix resp. rays"),
  ("vertex_indices", "Incidence Matrix resp. vertices"),
  ("vertex_and_ray_indices", "Incidence Matrix resp. vertices and rays"),
)
  M = Symbol(sym)
  _M = Symbol("_", sym)
  @eval begin
    $M(iter::SubObjectIterator) = $_M(Val(iter.Acc), iter.Obj; iter.options...)
    $_M(::Any, ::PolyhedralObjectUnion) =
      throw(ArgumentError(string($name, " not defined in this context.")))
  end
end

# Matrices with rational or integer elements
for (sym, name) in (
  ("point_matrix", "Point Matrix"),
  ("vector_matrix", "Vector Matrix"),
  ("generator_matrix", "Generator Matrix"),
)
  M = Symbol(sym)
  _M = Symbol("_", sym)
  @eval begin
    $M(iter::SubObjectIterator{<:AbstractVector{QQFieldElem}}) =
      matrix(QQ, $_M(Val(iter.Acc), iter.Obj; iter.options...))
    $M(iter::SubObjectIterator{<:AbstractVector{ZZRingElem}}) =
      matrix(ZZ, $_M(Val(iter.Acc), iter.Obj; iter.options...))
    $M(iter::SubObjectIterator{<:AbstractVector{<:FieldElem}}) =
      matrix(coefficient_field(iter.Obj), $_M(Val(iter.Acc), iter.Obj; iter.options...))
    $_M(::Any, ::PolyhedralObjectUnion) =
      throw(ArgumentError(string($name, " not defined in this context.")))
  end
end

function matrix_for_polymake(iter::SubObjectIterator; homogenized=false)
  if hasmethod(_matrix_for_polymake, Tuple{Val{iter.Acc}})
    return _matrix_for_polymake(Val(iter.Acc))(
      Val(iter.Acc), iter.Obj; homogenized=homogenized, iter.options...
    )
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

incidence_matrix(iter::SubObjectIterator) = IncidenceMatrix(iter)

# primitive generators only for ray based iterators
matrix(R::ZZRing, iter::SubObjectIterator{RayVector{QQFieldElem}}) =
  matrix(R, Polymake.common.primitive(matrix_for_polymake(iter)))
matrix(
  R::ZZRing, iter::SubObjectIterator{<:Union{RayVector{ZZRingElem},PointVector{ZZRingElem}}}
) = matrix(R, matrix_for_polymake(iter))
matrix(
  R::QQField,
  iter::SubObjectIterator{<:Union{RayVector{QQFieldElem},PointVector{QQFieldElem}}},
) = matrix(R, matrix_for_polymake(iter))
matrix(
  K::NCRing,
  iter::SubObjectIterator{<:Union{RayVector{<:FieldElem},PointVector{<:FieldElem}}},
) = matrix(K, matrix_for_polymake(iter))

function linear_matrix_for_polymake(iter::SubObjectIterator)
  if hasmethod(_linear_matrix_for_polymake, Tuple{Val{iter.Acc}})
    return _linear_matrix_for_polymake(Val(iter.Acc))(
      Val(iter.Acc), iter.Obj; iter.options...
    )
  elseif hasmethod(_affine_matrix_for_polymake, Tuple{Val{iter.Acc}})
    res = _affine_matrix_for_polymake(Val(iter.Acc))(
      Val(iter.Acc), iter.Obj; iter.options...
    )
    @req iszero(res[:, 1]) "Input not linear."
    return res[:, 2:end]
  end
  throw(ArgumentError("Linear Matrix for Polymake not defined in this context."))
end

function affine_matrix_for_polymake(iter::SubObjectIterator)
  if hasmethod(_affine_matrix_for_polymake, Tuple{Val{iter.Acc}})
    return _affine_matrix_for_polymake(Val(iter.Acc))(
      Val(iter.Acc), iter.Obj; iter.options...
    )
  elseif hasmethod(_linear_matrix_for_polymake, Tuple{Val{iter.Acc}})
    return homogenize(
      _linear_matrix_for_polymake(Val(iter.Acc))(Val(iter.Acc), iter.Obj; iter.options...),
      0,
    )
  end
  throw(ArgumentError("Affine Matrix for Polymake not defined in this context."))
end

Polymake.convert_to_pm_type(::Type{SubObjectIterator{RayVector{T}}}) where {T} =
  Polymake.Matrix{T}
Polymake.convert_to_pm_type(::Type{SubObjectIterator{PointVector{T}}}) where {T} =
  Polymake.Matrix{T}
Base.convert(::Type{<:Polymake.Matrix}, iter::SubObjectIterator) =
  assure_matrix_polymake(matrix_for_polymake(iter; homogenized=true))

function homogenized_matrix(
  field, x::SubObjectIterator{<:PointVector}, v::Union{Number,scalar_types_extended}=1
)
  @req isone(v) "PointVectors can only be (re-)homogenized with parameter 1, please convert to a matrix first"
  return matrix_for_polymake(x; homogenized=true)
end
function homogenized_matrix(
  field, x::SubObjectIterator{<:RayVector}, v::Union{Number,scalar_types_extended}=0
)
  @req iszero(v) "RayVectors can only be (re-)homogenized with parameter 0, please convert to a matrix first"
  return matrix_for_polymake(x; homogenized=true)
end

function homogenized_matrix(
  field, x::AbstractVector{<:PointVector}, v::Union{Number,scalar_types_extended}=1
)
  @req isone(v) "PointVectors can only be (re-)homogenized with parameter 1, please convert to a matrix first"
  return stack(field, [homogenize(coefficient_field(x[i]), x[i], v) for i in 1:length(x)])
end
function homogenized_matrix(
  field, x::AbstractVector{<:RayVector}, v::Union{Number,scalar_types_extended}=0
)
  @req iszero(v) "RayVectors can only be (re-)homogenized with parameter 0, please convert to a matrix first"
  return stack(field, [homogenize(coefficient_field(x[i]), x[i], v) for i in 1:length(x)])
end

homogenized_matrix(field, ::SubObjectIterator, v::Union{Number,scalar_types_extended}) =
  throw(ArgumentError("Content of SubObjectIterator not suitable for homogenized_matrix."))

unhomogenized_matrix(x::SubObjectIterator{<:RayVector}) = matrix_for_polymake(x)

unhomogenized_matrix(x::AbstractVector{<:PointVector}) =
  throw(ArgumentError("unhomogenized_matrix only meaningful for RayVectors"))

_ambient_dim(x::SubObjectIterator) = Polymake.polytope.ambient_dim(pm_object(x.Obj))

################################################################################

# Lineality often causes certain collections to be empty;
# the following definition allows to easily construct a working empty SOI

_empty_access() = nothing

function _empty_subobjectiterator(::Type{T}, Obj::PolyhedralObjectUnion) where {T}
  return SubObjectIterator{T}(Obj, _empty_access, 0, NamedTuple())
end

for f in ("_point_matrix", "_vector_matrix", "_generator_matrix")
  M = Symbol(f)
  @eval begin
    function $M(::Val{_empty_access}, P::PolyhedralObjectUnion; homogenized=false)
      typename = Polymake.bigobject_eltype(pm_object(P))
      T = if typename == "OscarNumber"
        Polymake.OscarNumber
      else
        _scalar_type_to_polymake(scalar_type_to_oscar[typename])
      end
      return Polymake.Matrix{T}(
        undef, 0, Polymake.polytope.ambient_dim(pm_object(P)) + homogenized
      )
    end
  end
end

for f in ("_facet_indices", "_ray_indices", "_vertex_indices", "_vertex_and_ray_indices")
  M = Symbol(f)
  @eval begin
    $M(::Val{_empty_access}, P::PolyhedralObjectUnion) =
      return Polymake.IncidenceMatrix(0, Polymake.polytope.ambient_dim(P))
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
const AbstractCollection = Dict{UnionAll,Union}([
  (PointVector, AnyVecOrMat),
  (RayVector, AnyVecOrMat),
  (
    LinearHalfspace,
    Union{AbstractVector{<:Halfspace},SubObjectIterator{<:Halfspace},AnyVecOrMat},
  ),
  (
    LinearHyperplane,
    Union{AbstractVector{<:Hyperplane},SubObjectIterator{<:Hyperplane},AnyVecOrMat},
  ),
  (
    AffineHalfspace,
    Union{
      AbstractVector{<:Halfspace},SubObjectIterator{<:Halfspace},Tuple{<:AnyVecOrMat,<:Any}
    },
  ),
  (
    AffineHyperplane,
    Union{
      AbstractVector{<:Hyperplane},
      SubObjectIterator{<:Hyperplane},
      Tuple{<:AnyVecOrMat,<:Any},
    },
  ),
])

affine_matrix_for_polymake(f::scalar_type_or_field, x::Union{Halfspace,Hyperplane}) =
  stack(f, augment(f, normal_vector(x), -negbias(x)))
affine_matrix_for_polymake(
  f::scalar_type_or_field, x::AbstractVector{<:Union{Halfspace,Hyperplane}}
) =
  stack(f, [affine_matrix_for_polymake(f, x[i]) for i in 1:length(x)])
linear_matrix_for_polymake(f::scalar_type_or_field, x::Union{Halfspace,Hyperplane}) =
  negbias(x) == 0 ? stack(f, normal_vector(x)) : throw(ArgumentError("Input not linear."))
linear_matrix_for_polymake(
  f::scalar_type_or_field, x::AbstractVector{<:Union{Halfspace,Hyperplane}}
) =
  stack(f, [linear_matrix_for_polymake(f, x[i]) for i in 1:length(x)])
