import Polymake: IncidenceMatrix

@doc raw"""
     IncidenceMatrix

A matrix with boolean entries. Each row corresponds to a fixed element of a collection of mathematical objects and the same holds for the columns and a second (possibly equal) collection. A `1` at entry `(i, j)` is interpreted as an incidence between object `i` of the first collection and object `j` of the second one.

# Examples
Note that the input and print of an `IncidenceMatrix` lists the non-zero indices for each row.
```jldoctest
julia> IM = IncidenceMatrix([[1,2,3],[4,5,6]])
2×6 IncidenceMatrix
[1, 2, 3]
[4, 5, 6]


julia> IM[1, 2]
true

julia> IM[2, 3]
false

julia> IM[:, 4]
2-element SparseVectorBool
[2]
```
"""
IncidenceMatrix

number_of_rows(i::IncidenceMatrix) = Polymake.nrows(i)
number_of_columns(i::IncidenceMatrix) = Polymake.ncols(i)

number_of_rows(A::Polymake.Matrix) = Polymake.nrows(A)
number_of_columns(A::Polymake.Matrix) = Polymake.ncols(A)

@doc raw"""
     row(i::IncidenceMatrix, n::Int)

Return the indices where the `n`-th row of `i` is `true`, as a `Set{Int}`.

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2,3],[4,5,6]])
2×6 IncidenceMatrix
[1, 2, 3]
[4, 5, 6]


julia> row(IM, 2)
Set{Int64} with 3 elements:
  5
  4
  6
```
"""
row(i::IncidenceMatrix, n::Int) = convert(Set{Int}, Polymake.row(i, n))

@doc raw"""
     column(i::IncidenceMatrix, n::Int)

Return the indices where the `n`-th column of `i` is `true`, as a `Set{Int}`.

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2,3],[4,5,6]])
2×6 IncidenceMatrix
[1, 2, 3]
[4, 5, 6]


julia> column(IM, 5)
Set{Int64} with 1 element:
  2
```
"""
column(i::IncidenceMatrix, n::Int) = convert(Set{Int}, Polymake.col(i, n))

function assure_matrix_polymake(m::Union{AbstractMatrix{Any}, AbstractMatrix{FieldElem}})
    a, b = size(m)
    if a > 0
        i = findfirst(_cannot_convert_to_fmpq, m)
        t = i === nothing ? QQFieldElem : typeof(m[i])
        if t <: Union{Polymake.Rational, Polymake.QuadraticExtension{Polymake.Rational}, Polymake.OscarNumber, Float64}
            m = Polymake.Matrix{Polymake.convert_to_pm_type(t)}(m)
        else
            m = Polymake.Matrix{_scalar_type_to_polymake(t)}(m)
        end
    else
        m = Polymake.Matrix{Polymake.Rational}(undef, a, b)
    end
    return m
end

assure_matrix_polymake(m::AbstractMatrix{<:FieldElem}) = Polymake.Matrix{Polymake.OscarNumber}(m)

assure_matrix_polymake(m::MatElem) = Polymake.Matrix{_scalar_type_to_polymake(eltype(m))}(m)

assure_matrix_polymake(m::Union{Oscar.ZZMatrix, Oscar.QQMatrix, AbstractMatrix{<:Union{QQFieldElem, ZZRingElem, Base.Integer, Base.Rational, Polymake.Rational, Polymake.QuadraticExtension, Polymake.OscarNumber, Float64}}}) = m

assure_matrix_polymake(m::SubArray{T, 2, U, V, W}) where {T<:Union{Polymake.Rational, Polymake.QuadraticExtension, Float64}, U, V, W} = Polymake.Matrix{T}(m)

function assure_vector_polymake(v::Union{AbstractVector{Any}, AbstractVector{FieldElem}})
    i = findfirst(_cannot_convert_to_fmpq, v)
    T = i === nothing ? QQFieldElem : typeof(v[i])
    return Polymake.Vector{_scalar_type_to_polymake(T)}(v)
end

assure_vector_polymake(v::AbstractVector{<:FieldElem}) = Polymake.Vector{Polymake.OscarNumber}(v)

assure_vector_polymake(v::AbstractVector{<:Union{QQFieldElem, ZZRingElem, Base.Integer, Base.Rational, Polymake.Rational, Polymake.QuadraticExtension, Polymake.OscarNumber, Float64}}) = v

affine_matrix_for_polymake(x::Tuple{<:AnyVecOrMat, <:AbstractVector}) = augment(unhomogenized_matrix(x[1]), -Vector(assure_vector_polymake(x[2])))
affine_matrix_for_polymake(x::Tuple{<:AnyVecOrMat, <:Any}) = homogenized_matrix(x[1], -x[2])

_cannot_convert_to_fmpq(x::Any) = !hasmethod(convert, Tuple{Type{QQFieldElem}, typeof(x)})

linear_matrix_for_polymake(x::Union{Oscar.ZZMatrix, Oscar.QQMatrix, AbstractMatrix}) = assure_matrix_polymake(x)

linear_matrix_for_polymake(x::AbstractVector{<:AbstractVector}) = assure_matrix_polymake(stack(x...))

matrix_for_polymake(x::Union{Oscar.ZZMatrix, Oscar.QQMatrix, AbstractMatrix}) = assure_matrix_polymake(x)

number_of_rows(x::SubArray{T, 2, U, V, W}) where {T, U, V, W} = size(x, 1)

function Polymake.Matrix{Polymake.Rational}(x::Union{Oscar.QQMatrix,AbstractMatrix{Oscar.QQFieldElem}})
    res = Polymake.Matrix{Polymake.Rational}(size(x)...)
    for i in eachindex(x)
        res[i] = x[i]
    end
    return res
end

function Polymake.Matrix{Polymake.OscarNumber}(x::Union{MatElem, AbstractMatrix{<:FieldElem}})
    res = Polymake.Matrix{Polymake.OscarNumber}(size(x)...)
    for i in eachindex(x)
        res[i] = x[i]
    end
    return res
end

_isempty_halfspace(x::Pair{<:Union{Oscar.MatElem, AbstractMatrix}, Any}) = isempty(x[1])
_isempty_halfspace(x) = isempty(x)


Base.convert(::Type{Polymake.QuadraticExtension{Polymake.Rational}}, x::QQFieldElem) = Polymake.QuadraticExtension(convert(Polymake.Rational, x))

Base.convert(T::Type{<:Polymake.Matrix}, x::Union{ZZMatrix,QQMatrix}) = Base.convert(T, Matrix(x))

Base.convert(::Type{<:Polymake.Integer}, x::ZZRingElem) = GC.@preserve x return Polymake.new_integer_from_fmpz(x)

Base.convert(::Type{<:Polymake.Rational}, x::QQFieldElem) = GC.@preserve x return Polymake.new_rational_from_fmpq(x)

Base.convert(::Type{<:Polymake.Rational}, x::ZZRingElem) = GC.@preserve x return Polymake.new_rational_from_fmpz(x)

Base.convert(::Type{<:Polymake.Integer}, x::QQFieldElem) = GC.@preserve x return Polymake.new_integer_from_fmpq(x)

function Base.convert(::Type{ZZRingElem}, x::Polymake.Integer)
    res = ZZRingElem()
    GC.@preserve x Polymake.new_fmpz_from_integer(x, pointer_from_objref(res))
    return res
end

function Base.convert(::Type{QQFieldElem}, x::Polymake.Rational)
    res = QQFieldElem()
    GC.@preserve x Polymake.new_fmpq_from_rational(x, pointer_from_objref(res))
    return res
end

function Base.convert(::Type{QQFieldElem}, x::Polymake.Integer)
    res = QQFieldElem()
    GC.@preserve x Polymake.new_fmpq_from_integer(x, pointer_from_objref(res))
    return res
end

function Base.convert(::Type{ZZRingElem}, x::Polymake.Rational)
    res = ZZRingElem()
    GC.@preserve x Polymake.new_fmpz_from_rational(x, pointer_from_objref(res))
    return res
end

Base.convert(::Type{Polymake.OscarNumber}, x::FieldElem) = Polymake.OscarNumber(x)

(::Type{T})(x::Polymake.OscarNumber) where T<:FieldElem = convert(T, Polymake.unwrap(x))

Base.convert(::Type{Polymake.Matrix{Polymake.OscarNumber}}, x::MatElem{<:FieldElem}) = Polymake.Matrix{Polymake.OscarNumber}(x)

(R::QQField)(x::Polymake.Rational) = convert(QQFieldElem, x)
(Z::ZZRing)(x::Polymake.Rational) = convert(ZZRingElem, x)

function (NF::Hecke.EmbeddedNumField)(x::Polymake.QuadraticExtension{Polymake.Rational})
    g = Polymake.generating_field_elements(x)
    if g.r == 0 || g.b == 0
        return NF(convert(QQFieldElem, g.a))
    end
    isq = Hecke.is_quadratic_type(number_field(NF))
    @req isq[1] "Can not construct non-trivial QuadraticExtension in non-quadratic number field."
    @req isq[2] == base_field(number_field(NF))(g.r) "Source and target fields do not match."
    a = NF(basis(number_field(NF))[2])
    return convert(QQFieldElem, g.a) + convert(QQFieldElem, g.b) * a
end

(F::Field)(x::Polymake.Rational) = F(QQ(x))
(F::Field)(x::Polymake.OscarNumber) = F(Polymake.unwrap(x))

Polymake.convert_to_pm_type(::Type{Oscar.ZZMatrix}) = Polymake.Matrix{Polymake.Integer}
Polymake.convert_to_pm_type(::Type{Oscar.QQMatrix}) = Polymake.Matrix{Polymake.Rational}
Polymake.convert_to_pm_type(::Type{Oscar.ZZRingElem}) = Polymake.Integer
Polymake.convert_to_pm_type(::Type{Oscar.QQFieldElem}) = Polymake.Rational
Polymake.convert_to_pm_type(::Type{T}) where T<:FieldElem = Polymake.OscarNumber
Polymake.convert_to_pm_type(::Type{<:Oscar.MatElem}) = Polymake.Matrix{Polymake.OscarNumber}

function remove_zero_rows(A::AbstractMatrix)
    A[findall(x->!iszero(x),collect(eachrow(A))),:]
end
function remove_zero_rows(A::Oscar.MatElem)
    remove_zero_rows(Matrix(A))
end

# function remove_redundant_rows(A::Union{Oscar.MatElem,AbstractMatrix})
#     rindices = Polymake.Set{Polymake.to_cxx_type(Int64)}(1:size(A, 1))
#     for i in rindices
#         for j in rindices
#             if i == j
#                 continue
#             end
#             if A[i, :] == A[j, :]
#                 delete!(rindices, j)
#             end
#         end
#     end
#     return A[rindices]
# end

function augment(vec::AbstractVector, val)
    s = size(vec)
    res = similar(vec, (s[1] + 1,))
    res[1] = val
    res[2:end] = vec
    return assure_vector_polymake(res)
end

function augment(mat::AbstractMatrix, vec::AbstractVector)
    s = size(mat)
    res = similar(mat, promote_type(eltype(mat), eltype(vec)),(s[1], s[2] + 1))
    res[:, 1] = vec
    res[:, 2:end] = mat
    return assure_matrix_polymake(res)
end

homogenize(vec::AbstractVector, val::Number = 0) = augment(vec, val)
homogenize(mat::AbstractMatrix, val::Number = 1) = augment(mat, fill(val, size(mat, 1)))
homogenize(mat::MatElem, val::Number = 1) = homogenize(Matrix(mat), val)
homogenize(nothing,val::Number)=nothing
homogenized_matrix(x::Union{AbstractVecOrMat,MatElem,Nothing}, val::Number) = homogenize(x, val)
homogenized_matrix(x::AbstractVector, val::Number) = permutedims(homogenize(x, val))
homogenized_matrix(x::AbstractVector{<:AbstractVector}, val::Number) = stack((homogenize(x[i], val) for i in 1:length(x))...)

dehomogenize(vec::AbstractVector) = vec[2:end]
dehomogenize(mat::AbstractMatrix) = mat[:, 2:end]

unhomogenized_matrix(x::AbstractVector) = assure_matrix_polymake(stack(x))
unhomogenized_matrix(x::AbstractMatrix) = assure_matrix_polymake(x)
unhomogenized_matrix(x::MatElem) = Matrix(assure_matrix_polymake(x))
unhomogenized_matrix(x::AbstractVector{<:AbstractVector}) = unhomogenized_matrix(stack(x...))

"""
    stack(A::AbstractVecOrMat, B::AbstractVecOrMat)

Stacks `A` and `B` vertically. The difference to `vcat`is that `AbstractVector`s are always
interpreted as row vectors. Empty vectors are ignored.

# Examples

```
julia> stack([1, 2], [0, 0])
2×2 Matrix{Int64}:
 1  2
 0  0

julia> stack([1 2], [0 0])
2×2 Matrix{Int64}:
 1  2
 0  0

julia> stack([1 2], [0, 0])
2×2 Matrix{Int64}:
 1  2
 0  0

julia> stack([1, 2], [0 0])
2×2 Matrix{Int64}:
 1  2
 0  0

julia> stack([1 2], [])
1×2 Matrix{Int64}:
 1  2
```
"""
stack(A::AbstractMatrix, ::Nothing) = A
stack(::Nothing, B::AbstractMatrix) = B
stack(A::AbstractMatrix, B::AbstractMatrix) = [A; B]
stack(A::AbstractMatrix, B::AbstractVector) = isempty(B) ? A :  [A; permutedims(B)]
stack(A::AbstractVector, B::AbstractMatrix) = isempty(A) ? B : [permutedims(A); B]
stack(A::AbstractVector, B::AbstractVector) = isempty(A) ? B : [permutedims(A); permutedims(B)]
stack(A::AbstractVector, ::Nothing) = permutedims(A)
stack(::Nothing, B::AbstractVector) = permutedims(B)
stack(x, y, z...) = stack(stack(x, y), z...)
stack(x) = stack(x, nothing)
# stack(x::Union{QQMatrix, ZZMatrix}, ::Nothing) = x
#=
function stack(A::Vector{Polymake.Vector{Polymake.Rational}})
    if length(A)==2
        return stack(A[1],A[2])
    end
    M=stack(A[1],A[2])
    for i in 3:length(A)
        M=stack(M,A[i])
    end
    return M
end
=#

_ambient_dim(x::AbstractVector) = length(x)
_ambient_dim(x::AbstractMatrix) = size(x, 2)
_ambient_dim(x::AbstractVector{<:AbstractVector}) = _ambient_dim(x[1])
_ambient_dim(x::MatElem) = ncols(x)

"""
    decompose_vdata(A::AbstractMatrix)

Given a (homogeneous) polymake matrix split into vertices and rays and dehomogenize.
"""
function decompose_vdata(A::AbstractMatrix)
    vertex_indices = findall(!iszero, view(A, :, 1))
    ray_indices = findall(iszero, view(A, :, 1))
    return (vertices = A[vertex_indices, 2:end], rays = A[ray_indices, 2:end])
end

function decompose_hdata(A)
    (A = -A[:, 2:end], b = A[:, 1])
end

# TODO: different printing within oscar? if yes, implement the following method
# Base.show(io::IO, ::MIME"text/plain", I::IncidenceMatrix) = show(io, "text/plain", Matrix{Bool}(I))

####################################################################
# Prepare interface for objects to be defined later
####################################################################

abstract type PolyhedralObject{T} end

# several toric types need to inherit from abstract schemes and thus cannot
# inherit from PolyhedralObject, but they still need to behave like polyhedral objects
# for some operations.
const PolyhedralObjectUnion = Union{PolyhedralObject, NormalToricVarietyType}

@doc raw"""
    coefficient_field(P::Union{Polyhedron{T}, Cone{T}, PolyhedralFan{T}, PolyhedralComplex{T}) where T<:scalar_types

Return the parent `Field` of the coefficients of `P`.

# Examples
```jldoctest
julia> c = cross_polytope(2)
Polyhedron in ambient dimension 2

julia> coefficient_field(c)
Rational field
```
"""
coefficient_field(x::PolyhedralObject) =  x.parent_field
coefficient_field(x::PolyhedralObject{QQFieldElem}) = QQ

_get_scalar_type(::PolyhedralObject{T}) where T = T
_get_scalar_type(::NormalToricVarietyType) = QQFieldElem

################################################################################
######## Scalar types
################################################################################
const scalar_types = Union{FieldElem, Float64}

const scalar_type_to_oscar = Dict{String, Type}([("Rational", QQFieldElem),
                                                 ("QuadraticExtension<Rational>", Hecke.EmbeddedNumFieldElem{nf_elem}),
                                                 ("QuadraticExtension", Hecke.EmbeddedNumFieldElem{nf_elem}),
                                                 ("Float", Float64)])

const scalar_types_extended = Union{scalar_types, ZZRingElem}

const scalar_type_or_field = Union{Type{<:scalar_types}, Field}

_scalar_type_to_polymake(::Type{QQFieldElem}) = Polymake.Rational
_scalar_type_to_polymake(::Type{<:FieldElem}) = Polymake.OscarNumber
_scalar_type_to_polymake(::Type{Float64}) = Float64

####################################################################
# Parent Fields
####################################################################

function _embedded_quadratic_field(r::ZZRingElem)
    if iszero(r)
        R, = rationals_as_number_field()
        return Hecke.embedded_field(R, real_embeddings(R)[])
    else
        R, = quadratic_field(r)
        return Hecke.embedded_field(R, real_embeddings(R)[2])
    end
end

function _find_parent_field(::Type{T}, x, y...) where T <: scalar_types
    f = _find_parent_field(T, x)
    elem_type(f) == T && return f
    return _find_parent_field(T, y...)
end
function _find_parent_field(::Type{T}, x::AbstractArray{<:FieldElem}) where T <: scalar_types
    for el in x
        el isa T && return parent(el)
    end
    return QQ
end
function _find_parent_field(::Type{T}, x::MatElem{<:FieldElem}) where T <: scalar_types
    f = base_ring(x)
    elem_type(f) == T && return f
    return QQ
end
_find_parent_field(::Type{T}, x::Tuple{<:AnyVecOrMat, <:Any}) where T <: scalar_types = _find_parent_field(T, x...)
_find_parent_field(::Type{T}, x::FieldElem) where T <: scalar_types = x isa T ? parent(x) : QQ
_find_parent_field(::Type{T}, x::Number) where T <: scalar_types = QQ
_find_parent_field(::Type{T}) where T <: scalar_types = QQ
# _find_parent_field() = QQ
_find_parent_field(::Type{T}, x::AbstractArray{<:AbstractArray}) where T <: scalar_types = _find_parent_field(T, x...)
function _find_parent_field(::Type{T}, x::AbstractArray) where T <: scalar_types
    for el in x
        el isa T && return parent(el)
    end
    return QQ
end

_determine_parent_and_scalar(f::Union{Field, ZZRing}, x...) = (f, elem_type(f))
# isempty(x) => standard/trivial field?
function _determine_parent_and_scalar(::Type{T}, x...) where T <: scalar_types
    if T == QQFieldElem
        f = QQ
    elseif T == Float64
        f = AbstractAlgebra.Floats{Float64}()
    else
        pf = _find_parent_field(T, x...)
        f = pf == QQ ? throw(ArgumentError("Scalars of type $T require specification of a parent field. Please pass the desired Field instead of the type or have a $T contained in your input data.")) : pf
    end
    return (f, T)
end

function _detect_default_field(::Type{Hecke.EmbeddedNumFieldElem{nf_elem}}, p::Polymake.BigObject)
    # we only want to check existing properties
    f = x -> Polymake.exists(p, string(x))
    propnames = intersect(propertynames(p), [:INPUT_RAYS, :POINTS, :RAYS, :VERTICES, :VECTORS, :INPUT_LINEALITY, :LINEALITY_SPACE, :FACETS, :INEQUALITIES, :EQUATIONS, :LINEAR_SPAN, :AFFINE_HULL])
    i = findfirst(f, propnames)
    # find first QuadraticExtension with root != 0
    # or first OscarNumber wrapping an embedded number field element
    while !isnothing(i)
        prop = getproperty(p, propnames[i])
        if eltype(prop) <: Polymake.QuadraticExtension
            for el in prop
                r = Polymake.generating_field_elements(el).r
                iszero(r) || return _embedded_quadratic_field(ZZ(r))[1]
            end
        elseif eltype(prop) <: Polymake.OscarNumber
            for el in prop
                on = Polymake.unwrap(el)
                if on isa Hecke.EmbeddedNumFieldElem{nf_elem}
                    return parent(on)
                end
            end
        end
        i = findnext(f, propnames, i + 1)
    end
    return _embedded_quadratic_field(ZZ(0))[1]
end

_detect_default_field(::Type{QQFieldElem}, p::Polymake.BigObject) = QQ
_detect_default_field(::Type{Float64}, p::Polymake.BigObject) = AbstractAlgebra.Floats{Float64}()

function _detect_default_field(::Type{T}, p::Polymake.BigObject) where T<:FieldElem
  # we only want to check existing properties
  propnames = intersect(Polymake.list_properties(p), ["INPUT_RAYS", "POINTS", "RAYS", "VERTICES", "VECTORS", "INPUT_LINEALITY", "LINEALITY_SPACE", "FACETS", "INEQUALITIES", "EQUATIONS", "LINEAR_SPAN", "AFFINE_HULL"])
  # find first OscarNumber wrapping a FieldElem
  for pn in propnames
    prop = getproperty(p, convert(String, pn))
    for el in prop
      on = Polymake.unwrap(el)
      if on isa T
        return parent(on)
      end
    end
  end
  throw(ArgumentError("BigObject does not contain information about a parent Field"))
end

function _detect_wrapped_type_and_field(p::Polymake.BigObject)
  # we only want to check existing properties
  propnames = intersect(Polymake.list_properties(p), ["INPUT_RAYS", "POINTS", "RAYS", "VERTICES", "VECTORS", "INPUT_LINEALITY", "LINEALITY_SPACE", "FACETS", "INEQUALITIES", "EQUATIONS", "LINEAR_SPAN", "AFFINE_HULL"])
  # find first OscarNumber wrapping a FieldElem
  for pn in propnames
    prop = getproperty(p, convert(String, pn))
    for el in prop
      on = Polymake.unwrap(el)
      if on isa FieldElem
        f = parent(on)
        T = elem_type(f)
        return (T, f)
      end
    end
  end
  throw(ArgumentError("BigObject does not contain information about a parent Field"))
end

function _detect_scalar_and_field(::Type{U}, p::Polymake.BigObject) where U<:PolyhedralObject
  T = detect_scalar_type(U, p)
  if isnothing(T)
    return _detect_wrapped_type_and_field(p)
  else
    return (T, _detect_default_field(T, p))
  end
end

# promotion helpers
function _promote_scalar_field(f::Union{Field, ZZRing}...)
  try
    x = sum([g(0) for g in f])
    p = parent(x)
    return (elem_type(p), p)
  catch e
    throw(ArgumentError("Can not find a mutual parent field for $f."))
    end
end

function _promote_scalar_field(a::AbstractArray{<:FieldElem})
    isempty(a) && return (QQFieldElem, QQ)
    return _promote_scalar_field(parent.(a)...)
end

_parent_or_coefficient_field(r::Base.RefValue{<:Union{FieldElem, ZZRingElem}}) = parent(r.x)

_parent_or_coefficient_field(v::AbstractVector{T}) where T<:Union{FieldElem, ZZRingElem} = _determine_parent_and_scalar(T, v)[1]

function _promoted_bigobject(::Type{T}, obj::PolyhedralObject{U}) where {T <: scalar_types, U <: scalar_types}
  T == U ? pm_object(obj) : Polymake.common.convert_to{_scalar_type_to_polymake(T)}(pm_object(obj))
end

# oscarnumber helpers

function Polymake._fieldelem_to_rational(e::EmbeddedElem)
   return Rational{BigInt}(QQ(e))
end

function Polymake._fieldelem_is_rational(e::EmbeddedElem)
   return is_rational(e)
end

function Polymake._fieldelem_to_float(e::EmbeddedElem)
   return Float64(real(embedding(parent(e))(data(e), 32)))
end

# convert a Polymake.BigObject's scalar from QuadraticExtension to OscarNumber (Polytope only)

function _polyhedron_qe_to_on(x::Polymake.BigObject, f::Field)
  res = Polymake.polytope.Polytope{Polymake.OscarNumber}()
  for pn in Polymake.list_properties(x)
    prop = Polymake.give(x, pn)
    Polymake.take(res, string(pn), _property_qe_to_on(prop, f))
  end
  return res
end

_property_qe_to_on(x::Polymake.BigObject, f::Field) = Polymake.BigObject(Polymake.bigobject_type(x), x)

_property_qe_to_on(x::Polymake.PropertyValue, f::Field) = x

_property_qe_to_on(x::Polymake.QuadraticExtension{Polymake.Rational}, f::Field) = f(x)

function _property_qe_to_on(x, f::Field)
  if hasmethod(length, (typeof(x),)) && eltype(x) <: Polymake.QuadraticExtension{Polymake.Rational}
    return f.(x)
  else
    return x
  end
end
