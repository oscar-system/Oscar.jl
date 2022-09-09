import Polymake: IncidenceMatrix

const nf_scalar = Union{nf_elem, fmpq}

function assure_matrix_polymake(m::Union{AbstractMatrix{Any}, AbstractMatrix{FieldElem}})
    a, b = size(m)
    if a > 0
        i = findfirst(_cannot_convert_to_fmpq, m)
        t = typeof(m[i])
        if t <: Union{values(scalar_type_to_polymake)...}
            m = Polymake.Matrix{Polymake.convert_to_pm_type(t)}(m)
        else
            m = Polymake.Matrix{scalar_type_to_polymake[t]}(m)
        end
    else
        m = Polymake.Matrix{Polymake.Rational}(undef, a, b)
    end
    return m
end

assure_matrix_polymake(m::AbstractMatrix{nf_scalar}) = Polymake.Matrix{Polymake.QuadraticExtension{Polymake.Rational}}(m)

assure_matrix_polymake(m::Union{Oscar.fmpz_mat, Oscar.fmpq_mat, AbstractMatrix{<:Union{fmpq, fmpz, Base.Integer, Base.Rational, Polymake.Rational, Polymake.QuadraticExtension, Float64}}}) = m

function assure_vector_polymake(v::Union{AbstractVector{Any}, AbstractVector{FieldElem}})
    i = findfirst(_cannot_convert_to_fmpq, v)
    v = Polymake.Vector{scalar_type_to_polymake[typeof(v[i])]}(v)
    return v
end

assure_vector_polymake(v::AbstractVector{nf_scalar}) = Polymake.Vector{Polymake.QuadraticExtension{Polymake.Rational}}(v)

assure_vector_polymake(v::AbstractVector{<:Union{fmpq, fmpz, nf_elem, Base.Integer, Base.Rational, Polymake.Rational, Polymake.QuadraticExtension}}) = v

affine_matrix_for_polymake(x::Tuple{<:AnyVecOrMat, <:AbstractVector}) = augment(unhomogenized_matrix(x[1]), -Vector(assure_vector_polymake(x[2])))
affine_matrix_for_polymake(x::Tuple{<:AnyVecOrMat, <:Any}) = homogenized_matrix(x[1], -x[2])

_cannot_convert_to_fmpq(x::Any) = !hasmethod(convert, Tuple{Type{fmpq}, typeof(x)})

linear_matrix_for_polymake(x::Union{Oscar.fmpz_mat, Oscar.fmpq_mat, AbstractMatrix}) = assure_matrix_polymake(x)

matrix_for_polymake(x::Union{Oscar.fmpz_mat, Oscar.fmpq_mat, AbstractMatrix}) = assure_matrix_polymake(x)

function Polymake.Matrix{Polymake.Rational}(x::Union{Oscar.fmpq_mat,AbstractMatrix{Oscar.fmpq}})
    res = Polymake.Matrix{Polymake.Rational}(size(x)...)
    for i in eachindex(x)
        res[i] = x[i]
    end
    return res
end

_isempty_halfspace(x::Pair{<:Union{Oscar.MatElem, AbstractMatrix}, Any}) = isempty(x[1])
_isempty_halfspace(x) = isempty(x)

Base.convert(::Type{Polymake.Integer}, x::fmpz) = Polymake.Integer(BigInt(x))
Base.convert(::Type{Polymake.Rational}, x::fmpz) = Polymake.Rational(convert(Polymake.Integer, x), convert(Polymake.Integer, 1))
Base.convert(::Type{Polymake.Rational}, x::fmpq) = Polymake.Rational(convert(Polymake.Integer, numerator(x)), convert(Polymake.Integer, denominator(x)))

# export nf_scalar

Base.zero(::Type{nf_scalar}) = fmpq()
# Base.one(::Type{nf_scalar}) = fmpq(1)

Base.convert(::Type{nf_scalar}, x::Number) = convert(fmpq, x)
Base.convert(::Type{nf_scalar}, x::nf_elem) = x

nf_scalar(x::Union{Number, nf_elem}) = convert(nf_scalar, x)

function Base.convert(::Type{Polymake.QuadraticExtension{Polymake.Rational}}, x::nf_elem)
    isq = Hecke.is_quadratic_type(parent(x))
    if !isq[1] || isq[2] < 0
        throw(ArgumentError("Conversion from nf_elem to QuadraticExtension{Rational} only defined for elements of real quadratic number fields defined by a polynomial of the form 'ax^2 - b'."))
    end
    r = convert(Polymake.Rational, isq[2])
    c = coordinates(x)
    return Polymake.QuadraticExtension{Polymake.Rational}(convert(Polymake.Rational, c[1]), convert(Polymake.Rational, c[2]), r)
end

Base.convert(::Type{Polymake.QuadraticExtension{Polymake.Rational}}, x::fmpq) = Polymake.QuadraticExtension(convert(Polymake.Rational, x))

function Base.convert(::Type{nf_scalar}, x::Polymake.QuadraticExtension{Polymake.Rational})
    g = Polymake.generating_field_elements(x)
    if g.r == 0 || g.b == 0
        return convert(fmpq, g.a)
    end
    R, a = quadratic_field(convert(fmpz,  g.r))
    return convert(fmpq, g.a) + convert(fmpq, g.b) * a
end

Base.convert(T::Type{<:Polymake.Matrix}, x::Union{fmpz_mat,fmpq_mat}) = Base.convert(T, Matrix(x))

Polymake.convert_to_pm_type(::Type{Oscar.fmpz_mat}) = Polymake.Matrix{Polymake.Integer}
Polymake.convert_to_pm_type(::Type{Oscar.fmpq_mat}) = Polymake.Matrix{Polymake.Rational}
Polymake.convert_to_pm_type(::Type{Oscar.fmpz}) = Polymake.Integer
Polymake.convert_to_pm_type(::Type{Oscar.fmpq}) = Polymake.Rational

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
# stack(x::Union{fmpq_mat, fmpz_mat}, ::Nothing) = x
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

Base.convert(::Type{fmpq}, q::Polymake.Rational) = fmpq(Polymake.numerator(q), Polymake.denominator(q))

Base.convert(::Type{fmpz}, q::Polymake.Rational) = convert(fmpz, convert(Polymake.Integer, q))

# TODO: different printing within oscar? if yes, implement the following method
# Base.show(io::IO, ::MIME"text/plain", I::IncidenceMatrix) = show(io, "text/plain", Matrix{Bool}(I))
