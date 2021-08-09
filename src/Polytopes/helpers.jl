matrix_for_polymake(x::Union{Oscar.fmpz_mat,AbstractMatrix{Oscar.fmpz}}) = Matrix{BigInt}(x)
matrix_for_polymake(x::Union{Oscar.fmpq_mat,AbstractMatrix{Oscar.fmpq}}) =
    Matrix{Hecke.Rational{BigInt}}(x)
matrix_for_polymake(x) = x

function remove_zero_rows(A::Union{Oscar.MatElem,AbstractMatrix})
    A[findall(x->!iszero(x),collect(eachrow(A))),:]
end

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


# This is a specific polymake data structure supporting fast functions
#  for rows->sets, rows containing col_i==true, etc.
struct IncidenceMatrix
   pm_incidencematrix::Polymake.IncidenceMatrix
end

function IncidenceMatrix(TrueIndices::Vector{Vector{Int64}})
   nrows = length(TrueIndices)
   ncols = maximum([maximum(set) for set in TrueIndices])
   IM = Polymake.IncidenceMatrix(nrows, ncols)
   i = 1
   for set in TrueIndices
      for j in set
         IM[i,j] = 1
      end
      i = i+1
  end
   return IncidenceMatrix(IM)
end



#TODO: change how incidence matrices are shown (not zero base but maybe bool?)
function Base.show(io::IO, I::IncidenceMatrix)
    show(io,"text/plain", (Matrix{Bool}(I.pm_incidencematrix)))
end
