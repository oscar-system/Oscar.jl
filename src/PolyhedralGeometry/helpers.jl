import Polymake: IncidenceMatrix

@doc raw"""
    IncidenceMatrix

A matrix with boolean entries. Each row corresponds to a fixed element of a collection of mathematical objects and the same holds for the columns and a second (possibly equal) collection. A `1` at entry `(i, j)` is interpreted as an incidence between object `i` of the first collection and object `j` of the second one.

# Examples
Note that the input of this example and the print of an `IncidenceMatrix` list the non-zero indices for each row.
```jldoctest
julia> IM = incidence_matrix([[1,2,3],[4,5,6]])
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

@doc raw"""
    incidence_matrix(r::Base.Integer, c::Base.Integer)

Return an `IncidenceMatrix` of size r x c whose entries are all `false`.

# Examples
```jldoctest
julia> IM = incidence_matrix(8, 5)
8×5 IncidenceMatrix
[]
[]
[]
[]
[]
[]
[]
[]

```
"""
incidence_matrix(r::Base.Integer, c::Base.Integer) = IncidenceMatrix(undef, r, c)

@doc raw"""
    incidence_matrix(mat::Union{AbstractMatrix{Bool}, IncidenceMatrix})

Convert `mat` to an `IncidenceMatrix`.

# Examples
```jldoctest
julia> IM = incidence_matrix([true false true false true false; false true false true false true])
2×6 IncidenceMatrix
[1, 3, 5]
[2, 4, 6]

```
"""
incidence_matrix(mat::Union{AbstractMatrix{Bool},IncidenceMatrix}) = IncidenceMatrix(mat)

@doc raw"""
    incidence_matrix(mat::AbstractMatrix)

Convert the `0`/`1` matrix `mat` to an `IncidenceMatrix`. Entries become `true` if the initial entry is `1` and `false` if the initial entry is `0`.

# Examples
```jldoctest
julia> IM = incidence_matrix([1 0 1 0 1 0; 0 1 0 1 0 1])
2×6 IncidenceMatrix
[1, 3, 5]
[2, 4, 6]

```
"""
function incidence_matrix(mat::AbstractMatrix)
  m, n = size(mat)
  for i in 1:m
    for j in 1:n
      iszero(mat[i, j]) || isone(mat[i, j]) ||
        throw(
          ArgumentError("incidence_matrix requires matrices with 0/1 or boolean entries.")
        )
    end
  end
  return IncidenceMatrix(mat)
end

@doc raw"""
    incidence_matrix(r::Base.Integer, c::Base.Integer, incidenceRows::AbstractVector{<:AbstractVector{<:Base.Integer}})

Return an `IncidenceMatrix` of size r x c. The i-th element of `incidenceRows` lists the indices of the `true` entries of the i-th row.

# Examples
```jldoctest
julia> IM = incidence_matrix(3, 4, [[2, 3], [1]])
3×4 IncidenceMatrix
[2, 3]
[1]
[]

```
"""
incidence_matrix(
  r::Base.Integer,
  c::Base.Integer,
  incidenceRows::AbstractVector{<:AbstractVector{<:Base.Integer}},
) = IncidenceMatrix(r, c, incidenceRows)

@doc raw"""
    incidence_matrix(incidenceRows::AbstractVector{<:AbstractVector{<:Base.Integer}})

Return an `IncidenceMatrix` where the i-th element of `incidenceRows` lists the indices of the `true` entries of the i-th row. The dimensions of the result are the smallest possible row and column count that can be deduced from the input.

# Examples
```jldoctest
julia> IM = incidence_matrix([[2, 3], [1]])
2×3 IncidenceMatrix
[2, 3]
[1]

```
"""
incidence_matrix(incidenceRows::AbstractVector{<:AbstractVector{<:Base.Integer}}) =
  IncidenceMatrix(incidenceRows)

number_of_rows(i::IncidenceMatrix) = Polymake.nrows(i)
number_of_columns(i::IncidenceMatrix) = Polymake.ncols(i)

number_of_rows(A::Polymake.Matrix) = Polymake.nrows(A)
number_of_columns(A::Polymake.Matrix) = Polymake.ncols(A)

@doc raw"""
    row(i::IncidenceMatrix, n::Int)

Return the indices where the `n`-th row of `i` is `true`, as a `Set{Int}`.

# Examples
```jldoctest
julia> IM = incidence_matrix([[1,2,3],[4,5,6]])
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
julia> IM = incidence_matrix([[1,2,3],[4,5,6]])
2×6 IncidenceMatrix
[1, 2, 3]
[4, 5, 6]


julia> column(IM, 5)
Set{Int64} with 1 element:
  2
```
"""
column(i::IncidenceMatrix, n::Int) = convert(Set{Int}, Polymake.col(i, n))

const _polymake_scalars = Union{
  Polymake.Integer,
  Polymake.Rational,
  Polymake.QuadraticExtension,
  Polymake.OscarNumber,
  Float64,
  Polymake.TropicalNumber,
}
const _polymake_compatible_scalars = Union{
  QQFieldElem,ZZRingElem,Base.Integer,Base.Rational,_polymake_scalars
}

function assure_matrix_polymake(m::Union{AbstractMatrix{Any},AbstractMatrix{FieldElem}})
  a, b = size(m)
  if a > 0
    i = findfirst(_cannot_convert_to_fmpq, m)
    t = i === nothing ? QQFieldElem : typeof(m[i])
    if t <: _polymake_scalars
      m = Polymake.Matrix{Polymake.convert_to_pm_type(t)}(m)
    else
      m = Polymake.Matrix{_scalar_type_to_polymake(t)}(m)
    end
  else
    m = Polymake.Matrix{Polymake.Rational}(undef, a, b)
  end
  return m
end

assure_matrix_polymake(m::AbstractMatrix{<:FieldElem}) =
  Polymake.Matrix{Polymake.OscarNumber}(m)

assure_matrix_polymake(m::MatElem) = Polymake.Matrix{_scalar_type_to_polymake(eltype(m))}(m)

assure_matrix_polymake(
  m::Union{Oscar.ZZMatrix,Oscar.QQMatrix,AbstractMatrix{<:_polymake_compatible_scalars}}
) = m

assure_matrix_polymake(m::SubArray{T,2,U,V,W}) where {T<:Union{_polymake_scalars},U,V,W} =
  Polymake.Matrix{T}(m)

function assure_vector_polymake(v::Union{AbstractVector{Any},AbstractVector{FieldElem}})
  i = findfirst(_cannot_convert_to_fmpq, v)
  T = i === nothing ? QQFieldElem : typeof(v[i])
  return Polymake.Vector{_scalar_type_to_polymake(T)}(v)
end

assure_vector_polymake(v::AbstractVector{<:FieldElem}) =
  Polymake.Vector{Polymake.OscarNumber}(v)

assure_vector_polymake(v::AbstractVector{<:_polymake_compatible_scalars}) = v

affine_matrix_for_polymake(x::Tuple{<:AnyVecOrMat,<:AbstractVector}) =
  augment(unhomogenized_matrix(x[1]), -Vector(assure_vector_polymake(x[2])))
affine_matrix_for_polymake(x::Tuple{<:AnyVecOrMat,<:Any}) = homogenized_matrix(x[1], -x[2])

_cannot_convert_to_fmpq(x::Any) = !hasmethod(convert, Tuple{Type{QQFieldElem},typeof(x)})

linear_matrix_for_polymake(x::Union{Oscar.ZZMatrix,Oscar.QQMatrix,AbstractMatrix}) =
  assure_matrix_polymake(x)

linear_matrix_for_polymake(x::AbstractVector{<:AbstractVector}) =
  assure_matrix_polymake(stack(x))

matrix_for_polymake(x::Union{Oscar.ZZMatrix,Oscar.QQMatrix,AbstractMatrix}) =
  assure_matrix_polymake(x)

number_of_rows(x::SubArray{T,2,U,V,W}) where {T,U,V,W} = size(x, 1)

_isempty_halfspace(x::Pair{<:Union{Oscar.MatElem,AbstractMatrix},Any}) = isempty(x[1])
_isempty_halfspace(x) = isempty(x)

function Polymake.Matrix{T}(
  x::Union{MatElem,AbstractMatrix{<:FieldElem}}
) where {T<:_polymake_scalars}
  res = Polymake.Matrix{T}(size(x)...)
  return res .= x
end

Base.convert(::Type{Polymake.Matrix{T}}, x::MatElem) where {T} = Polymake.Matrix{T}(x)

Base.convert(::Type{Polymake.QuadraticExtension{Polymake.Rational}}, x::QQFieldElem) =
  Polymake.QuadraticExtension(convert(Polymake.Rational, x))

Base.convert(::Type{<:Polymake.Integer}, x::ZZRingElem) =
  GC.@preserve x return Polymake.new_integer_from_fmpz(x)

Base.convert(::Type{<:Polymake.Rational}, x::QQFieldElem) =
  GC.@preserve x return Polymake.new_rational_from_fmpq(x)

Base.convert(::Type{<:Polymake.Rational}, x::ZZRingElem) =
  GC.@preserve x return Polymake.new_rational_from_fmpz(x)

Base.convert(::Type{<:Polymake.Integer}, x::QQFieldElem) =
  GC.@preserve x return Polymake.new_integer_from_fmpq(x)

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

(::Type{T})(x::Polymake.OscarNumber) where {T<:FieldElem} = convert(T, Polymake.unwrap(x))

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

# Disambiguation
(F::QQBarField)(x::Polymake.Rational) = F(QQ(x))
(F::QQBarField)(x::Polymake.OscarNumber) = F(Polymake.unwrap(x))

Polymake.convert_to_pm_type(::Type{typeof(min)}) = Polymake.Min
Polymake.convert_to_pm_type(::Type{typeof(max)}) = Polymake.Max

Polymake.convert_to_pm_type(::Type{ZZMatrix}) = Polymake.Matrix{Polymake.Integer}
Polymake.convert_to_pm_type(::Type{QQMatrix}) = Polymake.Matrix{Polymake.Rational}
Polymake.convert_to_pm_type(::Type{ZZRingElem}) = Polymake.Integer
Polymake.convert_to_pm_type(::Type{QQFieldElem}) = Polymake.Rational
Polymake.convert_to_pm_type(::Type{T}) where {T<:FieldElem} = Polymake.OscarNumber
Polymake.convert_to_pm_type(::Type{<:MatElem{T}}) where {T} =
  Polymake.Matrix{Polymake.convert_to_pm_type(T)}
Polymake.convert_to_pm_type(::Type{<:Graph{T}}) where {T<:Union{Directed,Undirected}} =
  Polymake.Graph{T}

Base.convert(
  ::Type{<:Polymake.Graph{T}}, g::Graph{T}
) where {T<:Union{Directed,Undirected}} = Oscar.pm_object(g)

function remove_zero_rows(A::AbstractMatrix)
  A[findall(!iszero, collect(eachrow(A))), :]
end
function remove_zero_rows(A::AbstractMatrix{Float64})
  A[
    findall(
      x -> !isapprox(x, zero(x); atol=Polymake._get_global_epsilon()), collect(eachrow(A))
    ),
    :,
  ]
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
  @req s[1] > 0 "cannot homogenize empty vector"
  res = similar(vec, (s[1] + 1,))
  res[1] = val + zero(first(vec))
  res[2:end] = vec
  return assure_vector_polymake(res)
end

function augment(mat::AbstractMatrix, vec::AbstractVector)
  s = size(mat)
  res = similar(mat, promote_type(eltype(mat), eltype(vec)), (s[1], s[2] + 1))
  res[:, 1] = vec
  res[:, 2:end] = mat
  return assure_matrix_polymake(res)
end

homogenize(vec::AbstractVector, val::Number=0) = augment(vec, val)
homogenize(mat::AbstractMatrix, val::Number=1) = augment(mat, fill(val, size(mat, 1)))
homogenize(mat::MatElem, val::Number=1) = homogenize(Matrix(mat), val)
homogenize(nothing, val::Number) = nothing
homogenized_matrix(x::Union{AbstractVecOrMat,MatElem,Nothing}, val::Number) =
  homogenize(x, val)
homogenized_matrix(x::AbstractVector, val::Number) = permutedims(homogenize(x, val))
homogenized_matrix(x::AbstractVector{<:AbstractVector}, val::Number) =
  stack([homogenize(x[i], val) for i in 1:length(x)])

dehomogenize(vm::AbstractVecOrMat) = Polymake.call_function(:polytope, :dehomogenize, vm)

unhomogenized_matrix(x::AbstractVector) = assure_matrix_polymake(stack(x))
unhomogenized_matrix(x::AbstractMatrix) = assure_matrix_polymake(x)
unhomogenized_matrix(x::MatElem) = Matrix(assure_matrix_polymake(x))
unhomogenized_matrix(x::AbstractVector{<:AbstractVector}) =
  unhomogenized_matrix(stack(x))

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
stack(A::AbstractMatrix, B::AbstractVector) = isempty(B) ? A : [A; permutedims(B)]
stack(A::AbstractVector, B::AbstractMatrix) = isempty(A) ? B : [permutedims(A); B]
stack(A::AbstractVector, B::AbstractVector) =
  isempty(A) ? B : [permutedims(A); permutedims(B)]
stack(A::AbstractVector, ::Nothing) = permutedims(A)
stack(::Nothing, B::AbstractVector) = permutedims(B)
function stack(VV::AbstractVector{<:AbstractVector})
  @req length(VV) > 0 "at least one vector required"
  if VERSION >= v"1.9"
    permutedims(Base.stack(VV))
  else
    permutedims(reshape(collect(Iterators.flatten(VV)), length(VV[1]), length(VV)))
  end
end
stack(x, y, z...) = reduce(stack, z; init=stack(x, y))
stack(x) = stack(x, nothing)

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
  return (vertices=A[vertex_indices, 2:end], rays=A[ray_indices, 2:end])
end

function decompose_hdata(A)
  (A=-A[:, 2:end], b=A[:, 1])
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
const PolyhedralObjectUnion = Union{PolyhedralObject,NormalToricVarietyType}

@doc raw"""
    coefficient_field(P::Union{Polyhedron{T}, Cone{T}, PolyhedralFan{T}, PolyhedralComplex{T}) where T<:scalar_types

Return the parent `Field` of the coefficients of `P`.

# Examples
```jldoctest
julia> c = cross_polytope(2)
Polytope in ambient dimension 2

julia> coefficient_field(c)
Rational field
```
"""
coefficient_field(x::PolyhedralObject) = x.parent_field
coefficient_field(x::PolyhedralObject{QQFieldElem}) = QQ

_get_scalar_type(::PolyhedralObject{T}) where {T} = T
_get_scalar_type(::NormalToricVarietyType) = QQFieldElem

################################################################################
######## Scalar types
################################################################################
const scalar_types = Union{FieldElem,Float64}

const scalar_type_to_oscar = Dict{String,Type}([
  ("Rational", QQFieldElem),
  ("QuadraticExtension<Rational>", Hecke.EmbeddedNumFieldElem{AbsSimpleNumFieldElem}),
  ("QuadraticExtension", Hecke.EmbeddedNumFieldElem{AbsSimpleNumFieldElem}),
  ("Float", Float64),
])

const scalar_types_extended = Union{scalar_types,ZZRingElem}

const scalar_type_or_field = Union{Type{<:scalar_types},Field}

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

function _check_field_polyhedral(::Type{T}) where {T}
  @req !(T <: NumFieldElem) "Number fields must be embedded, e.g. via `embedded_number_field`."
  @req hasmethod(isless, (T, T)) "Field must be ordered and have `isless` method."
end

_find_elem_type(x::Nothing) = Any
_find_elem_type(x::Any) = typeof(x)
_find_elem_type(x::Type) = x
_find_elem_type(x::Polymake.Rational) = QQFieldElem
_find_elem_type(x::Polymake.Integer) = ZZRingElem
_find_elem_type(x::AbstractArray) = reshape(_find_elem_type.(x), :)
_find_elem_type(x::Tuple) = reduce(vcat, _find_elem_type.(x))
_find_elem_type(x::AbstractArray{<:AbstractArray}) =
  reduce(vcat, _find_elem_type.(x); init=[])
_find_elem_type(x::MatElem) = [elem_type(base_ring(x))]

function _guess_fieldelem_type(x...)
  types = filter(!=(Any), _find_elem_type(x))
  T = QQFieldElem
  for t in types
    if t == Float64
      return Float64
    elseif promote_type(t, T) != T
      T = t
    end
  end
  return T
end

_parent_or_coefficient_field(::Type{Float64}, x...) = AbstractAlgebra.Floats{Float64}()
_parent_or_coefficient_field(::Type{ZZRingElem}, x...) = ZZ

_parent_or_coefficient_field(r::Base.RefValue{<:Union{FieldElem,ZZRingElem}}, x...) =
  parent(r.x)
_parent_or_coefficient_field(v::AbstractArray{T}) where {T<:Union{FieldElem,ZZRingElem}} =
  _parent_or_coefficient_field(T, v)

# QQ is done as a special case here to avoid ambiguities
function _parent_or_coefficient_field(::Type{T}, e::Any) where {T<:FieldElem}
  hasmethod(parent, (typeof(e),)) && elem_type(parent(e)) <: T ? parent(e) : missing
end

function _parent_or_coefficient_field(::Type{T}, c::MatElem) where {T<:FieldElem}
  elem_type(base_ring(c)) <: T ? base_ring(c) : missing
end

function _parent_or_coefficient_field(::Type{T}, c::AbstractArray) where {T<:FieldElem}
  first([collect(skipmissing(_parent_or_coefficient_field.(Ref(T), c))); missing])
end

function _parent_or_coefficient_field(::Type{T}, x, y...) where {T<:scalar_types}
  for c in (x, y...)
    p = _parent_or_coefficient_field(T, c)
    if p !== missing && elem_type(p) <: T
      return p
    end
  end
  missing
end

function _parent_or_coefficient_field(::Type{T}, c::Tuple) where {T<:FieldElem}
  return _parent_or_coefficient_field(T, c...)
end

function _determine_parent_and_scalar(f::Union{Field,ZZRing}, x...)
  _check_field_polyhedral(elem_type(f))
  return (f, elem_type(f))
end

function _determine_parent_and_scalar(::Type{T}, x...) where {T<:scalar_types}
  T == QQFieldElem && return (QQ, QQFieldElem)
  p = _parent_or_coefficient_field(T, x...)
  @req p !== missing "Scalars of type $T require specification of a parent field. Please pass the desired Field instead of the type or have a $T contained in your input data."
  _check_field_polyhedral(elem_type(p))
  return (p, elem_type(p))
end

function _detect_default_field(
  ::Type{Hecke.EmbeddedNumFieldElem{AbsSimpleNumFieldElem}}, p::Polymake.BigObject
)
  # we only want to check existing properties
  f = x -> Polymake.exists(p, string(x))
  propnames = intersect(
    propertynames(p),
    [
      :INPUT_RAYS,
      :POINTS,
      :RAYS,
      :VERTICES,
      :VECTORS,
      :INPUT_LINEALITY,
      :LINEALITY_SPACE,
      :FACETS,
      :INEQUALITIES,
      :EQUATIONS,
      :LINEAR_SPAN,
      :AFFINE_HULL,
    ],
  )
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
        if on isa Hecke.EmbeddedNumFieldElem{AbsSimpleNumFieldElem}
          return parent(on)
        end
      end
    end
    i = findnext(f, propnames, i + 1)
  end
  return _embedded_quadratic_field(ZZ(0))[1]
end

_detect_default_field(::Type{QQFieldElem}, p::Polymake.BigObject) = QQ
_detect_default_field(::Type{Float64}, p::Polymake.BigObject) =
  AbstractAlgebra.Floats{Float64}()

function _detect_default_field(::Type{T}, p::Polymake.BigObject) where {T<:FieldElem}
  # we only want to check existing properties
  propnames = intersect(
    Polymake.list_properties(p),
    [
      "INPUT_RAYS",
      "POINTS",
      "RAYS",
      "VERTICES",
      "VECTORS",
      "INPUT_LINEALITY",
      "LINEALITY_SPACE",
      "FACETS",
      "INEQUALITIES",
      "EQUATIONS",
      "LINEAR_SPAN",
      "AFFINE_HULL",
    ],
  )
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
  propnames = intersect(
    Polymake.list_properties(p),
    [
      "INPUT_RAYS",
      "POINTS",
      "RAYS",
      "VERTICES",
      "VECTORS",
      "INPUT_LINEALITY",
      "LINEALITY_SPACE",
      "FACETS",
      "INEQUALITIES",
      "EQUATIONS",
      "LINEAR_SPAN",
      "AFFINE_HULL",
    ],
  )
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

function _detect_scalar_and_field(
  ::Type{U}, p::Polymake.BigObject
) where {U<:PolyhedralObject}
  T = detect_scalar_type(U, p)
  if isnothing(T)
    return _detect_wrapped_type_and_field(p)
  else
    return (T, _detect_default_field(T, p))
  end
end

# promotion helpers
function _promote_scalar_field(f::Union{Field,ZZRing}...)
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

function _promoted_bigobject(
  ::Type{T}, obj::PolyhedralObject{U}
) where {T<:scalar_types,U<:scalar_types}
  if T == U
    pm_object(obj)
  else
    Polymake.common.convert_to{_scalar_type_to_polymake(T)}(pm_object(obj))
  end
end

# oscarnumber helpers

function Polymake._fieldelem_to_floor(e::Union{EmbeddedNumFieldElem,QQBarFieldElem})
  return BigInt(floor(ZZRingElem, e))
end

function Polymake._fieldelem_to_ceil(e::Union{EmbeddedNumFieldElem,QQBarFieldElem})
  return BigInt(ceil(ZZRingElem, e))
end

function Polymake._fieldelem_to_rational(e::EmbeddedNumFieldElem)
  return Rational{BigInt}(QQ(e))
end

function Polymake._fieldelem_is_rational(e::EmbeddedNumFieldElem)
  return is_rational(e)
end

function Polymake._fieldelem_to_float(e::EmbeddedNumFieldElem)
  return Float64(real(embedding(parent(e))(data(e), 32)))
end

function Polymake._fieldelem_to_float(e::QQBarFieldElem)
  return Float64(ArbField(64)(e))
end

function Polymake._fieldelem_from_rational(::QQBarField, r::Rational{BigInt})
  return QQBarFieldElem(QQFieldElem(r))
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

_property_qe_to_on(x::Polymake.BigObject, f::Field) =
  Polymake.BigObject(Polymake.bigobject_type(x), x)

_property_qe_to_on(x::Polymake.PropertyValue, f::Field) = x

_property_qe_to_on(x::Polymake.QuadraticExtension{Polymake.Rational}, f::Field) = f(x)

function _property_qe_to_on(x, f::Field)
  if hasmethod(length, (typeof(x),)) &&
    eltype(x) <: Polymake.QuadraticExtension{Polymake.Rational}
    return f.(x)
  else
    return x
  end
end

# Helper function for conversion
# Cone -> Polyhedron
# PolyhedralFan -> PolyhedralComplex
# for transforming coordinate matrices.
function embed_at_height_one(M::AbstractMatrix{T}, add_vert::Bool) where {T}
  result = Polymake.Matrix{Polymake.convert_to_pm_type(T)}(
    add_vert + size(M)[1], size(M)[2] + 1
  )
  if add_vert
    result[1, 1] = 1
  end
  result[(add_vert + 1):end, 2:end] = M
  return result
end
