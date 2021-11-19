import Polymake: IncidenceMatrix

import Base: convert


matrix_for_polymake(x::Union{Oscar.fmpz_mat,AbstractMatrix{Oscar.fmpz}}) = Matrix{BigInt}(x)
matrix_for_polymake(x::Union{Oscar.fmpq_mat,AbstractMatrix{Oscar.fmpq}}) =
    Polymake.Matrix{Polymake.Rational}(x)
matrix_for_polymake(x) = x

affine_matrix_for_polymake(x::Tuple) = matrix_for_polymake(hcat(-x[2], x[1]))

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

function remove_zero_rows(A::Union{Oscar.MatElem,AbstractMatrix})
    A[findall(x->!iszero(x),collect(eachrow(A))),:]
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
    return res
end

function augment(mat::AbstractMatrix, vec::AbstractVector)
    s = size(mat)
    res = similar(mat, (s[1], s[2] + 1))
    res[:, 1] = vec
    res[:, 2:end] = mat
    return res
end

homogenize(vec::AbstractVector, val::Number = 0) = augment(vec, val)
homogenize(mat::AbstractMatrix, val::Number = 1) = augment(mat, fill(val, size(mat, 1)))
homogenize(mat::MatElem, val::Number = 1) = homogenize(Matrix(mat), val)
homogenize(nothing,val::Number)=nothing

dehomogenize(vec::AbstractVector) = vec[2:end]
dehomogenize(mat::AbstractMatrix) = mat[:, 2:end]

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
stack(A::AbstractMatrix,nothing) = A
stack(nothing,B::AbstractMatrix) = B
stack(A::AbstractMatrix, B::AbstractMatrix) = [A; B]
stack(A::AbstractMatrix, B::AbstractVector) = isempty(B) ? A :  [A; B']
stack(A::AbstractVector, B::AbstractMatrix) = isempty(A) ? B : [A'; B]
stack(A::AbstractVector, B::AbstractVector) = isempty(A) ? B : [A'; B']
stack(A::AbstractVector,nothing) = A'
stack(nothing,B::AbstractVector) = B'
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

# TODO: different printing within oscar? if yes, implement the following method
# Base.show(io::IO, ::MIME"text/plain", I::IncidenceMatrix) = show(io, "text/plain", Matrix{Bool}(I))
