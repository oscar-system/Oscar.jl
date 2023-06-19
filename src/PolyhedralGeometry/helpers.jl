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

nrows(i::IncidenceMatrix) = Polymake.nrows(i)
ncols(i::IncidenceMatrix) = Polymake.ncols(i)

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
        t = typeof(m[i])
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

assure_matrix_polymake(m::MatElem) = Polymake.OscarNumber.(m)

assure_matrix_polymake(m::Union{Oscar.ZZMatrix, Oscar.QQMatrix, AbstractMatrix{<:Union{QQFieldElem, ZZRingElem, Base.Integer, Base.Rational, Polymake.Rational, Polymake.QuadraticExtension, Polymake.OscarNumber, Float64}}}) = m

assure_matrix_polymake(m::SubArray{T, 2, U, V, W}) where {T<:Union{Polymake.Rational, Polymake.QuadraticExtension, Float64}, U, V, W} = Polymake.Matrix{T}(m)

function assure_vector_polymake(v::Union{AbstractVector{Any}, AbstractVector{FieldElem}})
    i = findfirst(_cannot_convert_to_fmpq, v)
    v = Polymake.Vector{_scalar_type_to_polymake(typeof(v[i]))}(v)
    return v
end

assure_vector_polymake(v::AbstractVector{<:FieldElem}) = Polymake.Vector{Polymake.OscarNumber}(v)

assure_vector_polymake(v::AbstractVector{<:Union{QQFieldElem, ZZRingElem, nf_elem, Base.Integer, Base.Rational, Polymake.Rational, Polymake.QuadraticExtension, Polymake.OscarNumber, Float64}}) = v

affine_matrix_for_polymake(x::Tuple{<:AnyVecOrMat, <:AbstractVector}) = augment(unhomogenized_matrix(x[1]), -Vector(assure_vector_polymake(x[2])))
affine_matrix_for_polymake(x::Tuple{<:AnyVecOrMat, <:Any}) = homogenized_matrix(x[1], -x[2])

_cannot_convert_to_fmpq(x::Any) = !hasmethod(convert, Tuple{Type{QQFieldElem}, typeof(x)})

linear_matrix_for_polymake(x::Union{Oscar.ZZMatrix, Oscar.QQMatrix, AbstractMatrix}) = assure_matrix_polymake(x)

linear_matrix_for_polymake(x::AbstractVector{<:AbstractVector}) = assure_matrix_polymake(stack(x...))

matrix_for_polymake(x::Union{Oscar.ZZMatrix, Oscar.QQMatrix, AbstractMatrix}) = assure_matrix_polymake(x)

nrows(x::SubArray{T, 2, U, V, W}) where {T, U, V, W} = size(x, 1)

function Polymake.Matrix{Polymake.Rational}(x::Union{Oscar.QQMatrix,AbstractMatrix{Oscar.QQFieldElem}})
    res = Polymake.Matrix{Polymake.Rational}(size(x)...)
    for i in eachindex(x)
        res[i] = x[i]
    end
    return res
end

_isempty_halfspace(x::Pair{<:Union{Oscar.MatElem, AbstractMatrix}, Any}) = isempty(x[1])
_isempty_halfspace(x) = isempty(x)

function Base.convert(::Type{Polymake.QuadraticExtension{Polymake.Rational}}, x::nf_elem)
    isq = Hecke.is_quadratic_type(parent(x))
    @req isq[1] && isq[2] >= 0 "Conversion from nf_elem to QuadraticExtension{Rational} only defined for elements of real quadratic number fields defined by a polynomial of the form 'ax^2 - b'"
    r = convert(Polymake.Rational, isq[2])
    c = coordinates(x)
    return Polymake.QuadraticExtension{Polymake.Rational}(convert(Polymake.Rational, c[1]), convert(Polymake.Rational, c[2]), r)
end

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

Base.convert(T::Type{<:FieldElem}, x::Polymake.OscarNumber) = convert(T, Polymake.unwrap(x))

(R::QQField)(x::Polymake.Rational) = convert(QQFieldElem, x)

function (NF::AnticNumberField)(x::Polymake.QuadraticExtension{Polymake.Rational})
    g = Polymake.generating_field_elements(x)
    if g.r == 0 || g.b == 0
        return NF(convert(QQFieldElem, g.a))
    end
    isq = Hecke.is_quadratic_type(NF)
    @req isq[2] == base_field(NF)(g.r) "Source and target fields do not match."
    a = basis(NF)[2]
    return convert(QQFieldElem, g.a) + convert(QQFieldElem, g.b) * a
end

(F::Field)(x::Polymake.Rational) = F(QQ(x))
(F::Field)(x::Polymake.OscarNumber) = F(Polymake.unwrap(x))

Polymake.convert_to_pm_type(::Type{Oscar.ZZMatrix}) = Polymake.Matrix{Polymake.Integer}
Polymake.convert_to_pm_type(::Type{Oscar.QQMatrix}) = Polymake.Matrix{Polymake.Rational}
Polymake.convert_to_pm_type(::Type{Oscar.ZZRingElem}) = Polymake.Integer
Polymake.convert_to_pm_type(::Type{Oscar.QQFieldElem}) = Polymake.Rational
Polymake.convert_to_pm_type(::Type{T}) where T<:FieldElem = Polymake.OscarNumber
# Polymake.convert_to_pm_type(::Type{Oscar.MatElem}) = Polymake.Matrix{Polymake.OscarNumber}

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
# Parent Fields
####################################################################

_common_parent(x::Field, y::Field) = x == QQ ? y : x
_find_parent_field(x, y...) = _common_parent(_find_parent_field(x), _find_parent_field(y...))
_find_parent_field(x::AbstractArray{<:FieldElem}) = length(x) == 0 ? QQ : parent(iterate(x)[1])
_find_parent_field(x::MatElem{<:FieldElem}) = length(x) == 0 ? QQ : base_ring(x)
_find_parent_field(x::Tuple{<:AnyVecOrMat, <:Any}) = _find_parent_field(x...)
_find_parent_field(x::FieldElem) = parent(x)
_find_parent_field(x::Number) = QQ
_find_parent_field(::Nothing) = QQ
_find_parent_field() = QQ
# fallback
function _find_parent_field(x::AbstractArray)
    length(x) == 0 && return QQ
    fe = findall(e -> e isa FieldElem, x)
    for i in fe
        parent(x[i]) == QQ || return parent(x[i])
    end
    return QQ
end

_determine_parent_and_scalar(f::Field, x...) = (f, elem_type(f))
# isempty(x) => standard/trivial field?
function _determine_parent_and_scalar(::Type{T}, x...) where T <: FieldElem
    if T == QQFieldElem
        f = QQ
    else
        pf = _find_parent_field(x...)
        f = pf == QQ ? throw(ArgumentError("Scalars of type $T require specification of a parent field. Please pass the desired Field instead of the type or have a $T contained in your input data.")) : pf
    end
    return (f, T)
end
