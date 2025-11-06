
# TODO : in this file, some functions for matrices and vectors are defined just to make other files work,
# such as forms.jl, transform_form.jl, linear_conjugate.jl and linear_centralizer.jl
# TODO : functions in this file are only temporarily, and often inefficient.
# TODO: once similar working methods are defined in other files or packages (e.g. Hecke),
# functions in this file are to be removed / moved / replaced
# TODO: when this happens, files mentioned above need to be modified too.




########################################################################
#
# Matrix manipulation
#
########################################################################

@doc raw"""
    permutation_matrix(R::NCRing, Q::AbstractVector{T}) where T <: Int
    permutation_matrix(R::NCRing, p::PermGroupElem)

Return the permutation matrix over the ring `R` corresponding to the sequence `Q` or to the permutation `p`.
If `Q` is a sequence, then `Q` must contain exactly once every integer from 1 to some `n`.

# Examples
```jldoctest
julia> s = perm([3, 1, 2])
(1,3,2)

julia> permutation_matrix(QQ, s)
[0   0   1]
[1   0   0]
[0   1   0]
```
"""
function permutation_matrix(R::NCRing, Q::AbstractVector{<:IntegerUnion}; check::Bool = true)
  if check
    @assert Set(Q) == Set(1:length(Q)) "$Q must be a permutation"
  end
  Z = zero_matrix(R, length(Q), length(Q))
  for i in 1:length(Q)
    Z[i, Q[i]] = 1
  end
  return Z
end

permutation_matrix(R::NCRing, p::PermGroupElem) = permutation_matrix(R, Vector(p), check = false)


########################################################################
#
# New operations
#
########################################################################

Base.:*(v::AbstractAlgebra.Generic.FreeModuleElem{T},x::MatrixGroupElem{T}) where T <: RingElem = v.parent(v.v*matrix(x))

Base.:*(v::Vector{T}, x::MatrixGroupElem{T}) where T <: RingElem = v*matrix(x)
Base.:*(x::MatrixGroupElem{T}, u::Vector{T}) where T <: RingElem = matrix(x)*u

# `on_tuples` and `on_sets` delegate to an action via `^` on the subobjects
# (`^` is the natural action in GAP)
Base.:^(v::AbstractAlgebra.Generic.FreeModuleElem{T},x::MatrixGroupElem{T}) where T <: RingElem = v.parent(v.v*matrix(x))

# action of matrix group elements on subspaces of a vector space
function Base.:^(V::AbstractAlgebra.Generic.Submodule{T}, x::MatrixGroupElem{T}) where T <: RingElem
  return sub(V.m, [v^x for v in V.gens])[1]
end
