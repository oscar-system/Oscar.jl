
# TODO : in this file, some functions for matrices and vectors are defined just to make other files work,
# such as forms.jl, transform_form.jl, linear_conjugate.jl and linear_centralizer.jl
# TODO : functions in this file are only temporarily, and often inefficient.
# TODO: once similar working methods are defined in other files or packages (e.g. Hecke), 
# functions in this file are to be removed / moved / replaced
# TODO: when this happens, files mentioned above need to be modified too.

export
    complement,
    conjugate_transpose,
    is_hermitian_matrix,
    is_skewsymmetric_matrix,
    lower_triangular_matrix,
    permutation_matrix,
    upper_triangular_matrix



########################################################################
#
# Matrix manipulation
#
########################################################################


"""
    matrix(A::Vector{AbstractAlgebra.Generic.FreeModuleElem{T}})

Return the matrix whose rows are the vectors in `A`.
All vectors in `A` must have the same length and the same base ring.
"""
function matrix(A::Vector{AbstractAlgebra.Generic.FreeModuleElem{T}}) where T <: FieldElem
   c = length(A[1].v)
   @assert all(x -> length(x.v)==c, A) "Vectors must have the same length"
   X = zero_matrix(base_ring(A[1]), length(A), c)
   for i in 1:length(A), j in 1:c
      X[i,j] = A[i][j]
   end

   return X
end

@doc Markdown.doc"""
    upper_triangular_matrix(L)

Return the upper triangular matrix whose entries on and above the diagonal are the elements of `L`.

An exception is thrown whenever the length of `L` is not equal to $n(n+1)/2$,
for some integer $n$.
"""
function upper_triangular_matrix(L)
   T = eltype(L)
   @assert T <: RingElem "L must be a collection of ring elements"
   d = Int(floor((sqrt(1+8*length(L))-1)/2))
   length(L)==div(d*(d+1),2) || throw(ArgumentError("Input vector of invalid length"))
   R = parent(L[1])
   x = zero_matrix(R,d,d)
   pos=1
   for i in 1:d, j in i:d
      x[i,j] = L[pos]
      pos+=1
   end
   return x
end

@doc Markdown.doc"""
    lower_triangular_matrix(L)

Return the upper triangular matrix whose entries on and below the diagonal are the elements of `L`.

An exception is thrown whenever the length of `L` is not equal to $n(n+1)/2$,
for some integer $n$.
"""
function lower_triangular_matrix(L)
   T = eltype(L)
   @assert T <: RingElem "L must be a collection of ring elements"
   d = Int(floor((sqrt(1+8*length(L))-1)/2))
   length(L)==div(d*(d+1),2) || throw(ArgumentError("Input vector of invalid length"))
   R = parent(L[1])
   x = zero_matrix(R,d,d)
   pos=1
   for i in 1:d, j in 1:i
      x[i,j] = L[pos]
      pos+=1
   end
   return x
end

"""
    conjugate_transpose(x::MatElem{T}) where T <: FinFieldElem

If the base ring of `x` is `GF(q^2)`, return the matrix `transpose( map ( y -> y^q, x) )`.
 An exception is thrown if the base ring does not have even degree.
"""
function conjugate_transpose(x::MatElem{T}) where T <: FinFieldElem
   iseven(degree(base_ring(x))) || throw(ArgumentError("The base ring must have even degree"))
   e = div(degree(base_ring(x)),2)
   return transpose(map(y -> frobenius(y,e),x))
end


# computes a complement for W in V (i.e. a subspace U of V such that V is direct sum of U and W)
"""
    complement(V::AbstractAlgebra.Generic.FreeModule{T}, W::AbstractAlgebra.Generic.Submodule{T})

Return a complement for `W` in `V`, i.e. a subspace `U` of `V` such that `V` is direct sum of `U` and `W`.
"""
function complement(V::AbstractAlgebra.Generic.FreeModule{T}, W::AbstractAlgebra.Generic.Submodule{T}) where T <: FieldElem
   @assert is_submodule(V,W) "The second argument is not a subspace of the first one"
   if dim(W)==0 return sub(V,basis(V)) end

   e = W.map

   H = matrix( vcat([e(g) for g in gens(W)], [zero(V) for i in 1:(dim(V)-dim(W)) ]) )
   A_left = identity_matrix(base_ring(V), dim(V))
   A_right = identity_matrix(base_ring(V), dim(V))
   for rn in 1:dim(W)     # rn = row number
      cn = rn    # column number
      while H[rn,cn]==0 cn+=1 end   # bring on the left the first non-zero entry
      swap_cols!(H,rn,cn)
      swap_rows!(A_right,rn,cn)
      for j in rn+1:dim(W)
         add_row!(H,H[j,rn]*H[rn,rn]^-1,rn,j)
         add_column!(A_left,A_left[j,rn]*A_left[rn,rn]^-1,j,rn)
      end
   end
   for j in dim(W)+1:dim(V)  H[j,j]=1  end
   H = A_left*H*A_right
   _gens = [V([H[i,j] for j in 1:dim(V)]) for i in dim(W)+1:dim(V) ]

   return sub(V,_gens)
end

"""
    permutation_matrix(F::Ring, Q::AbstractVector{T}) where T <: Int
    permutation_matrix(F::Ring, p::PermGroupElem)

Return the permutation matrix over the ring `R` corresponding to the sequence `Q` or to the permutation `p`.
If `Q` is a sequence, then `Q` must contain exactly once every integer from 1 to some `n`.
"""
function permutation_matrix(F::Ring, Q::AbstractVector{<:IntegerUnion})
   @assert Set(Q)==Set(1:length(Q)) "Invalid input"
   Z = zero_matrix(F,length(Q),length(Q))
   for i in 1:length(Q) Z[i,Q[i]] = 1 end
   return Z
end

permutation_matrix(F::Ring, p::PermGroupElem) = permutation_matrix(F, Vector(p))

^(a::MatElem, b::fmpz) = Hecke._generic_power(a, b)

########################################################################
#
# New properties
#
########################################################################

# TODO: not sure whether this definition of skew-symmetric is standard (for fields of characteristic 2)
"""
    is_skewsymmetric_matrix(B::MatElem{T}) where T <: Ring

Return whether the matrix `B` is skew-symmetric,
i.e. `B = -transpose(B)` and `B` has zeros on the diagonal.
Return `false` if `B` is not a square matrix.
"""
function is_skewsymmetric_matrix(B::MatElem{T}) where T <: RingElem
   n = nrows(B)
   n==ncols(B) || return false

   for i in 1:n
      B[i,i]==0 || return false
      for j in i+1:n
         B[i,j]==-B[j,i] || return false
      end
   end

   return true
end

"""
    is_hermitian_matrix(B::MatElem{T}) where T <: FinFieldElem

Return whether the matrix `B` is hermitian, i.e. `B = conjugate_transpose(B)`.
Return `false` if `B` is not a square matrix, or the field has not even degree.
"""
function is_hermitian_matrix(B::MatElem{T}) where T <: FinFieldElem
   n = nrows(B)
   n==ncols(B) || return false
   e = degree(base_ring(B))
   iseven(e) ? e = div(e,2) : return false

   for i in 1:n
      for j in i:n
         B[i,j]==frobenius(B[j,i],e) || return false
      end
   end

   return true
end

# return (true, h) if y = hx, (false, nothing) otherwise
# FIXME: at the moment, works only for fields
function _is_scalar_multiple_mat(x::MatElem{T}, y::MatElem{T}) where T <: RingElem
   F=base_ring(x)
   F==base_ring(y) || return (false, nothing)
   nrows(x)==nrows(y) || return (false, nothing)
   ncols(x)==ncols(y) || return (false, nothing)

   for i in 1:nrows(x), j in 1:ncols(x)
      if !iszero(x[i,j])
         h = y[i,j] * x[i,j]^-1
         return y == h*x ? (true,h) : (false, nothing)
      end
   end
  
   # at this point, x must be zero
   return y == 0 ? (true, F(1)) : (false, nothing)
end

########################################################################
#
# New operations
#
########################################################################

Base.getindex(V::AbstractAlgebra.Generic.FreeModule, i::Int) = gen(V, i)

# scalar product
Base.:*(v::AbstractAlgebra.Generic.FreeModuleElem{T},u::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: FieldElem = (v.v*transpose(u.v))[1]

Base.:*(v::AbstractAlgebra.Generic.FreeModuleElem{T},x::MatElem{T}) where T <: FieldElem = v.parent(v.v*x)
Base.:*(x::MatElem{T},u::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: FieldElem = x*transpose(u.v)

# evaluation of the form x into the vectors v and u
Base.:*(v::AbstractAlgebra.Generic.FreeModuleElem{T},x::MatElem{T},u::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: FieldElem = (v.v*x*transpose(u.v))[1]


Base.:*(v::AbstractAlgebra.Generic.FreeModuleElem{T},x::MatrixGroupElem{T}) where T <: FieldElem = v.parent(v.v*x.elm)
Base.:*(x::MatrixGroupElem{T},u::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: FieldElem = x.elm*transpose(u.v)

# `on_tuples` and `on_sets` delegate to an action via `^` on the subobjects
# (`^` is the natural action in GAP)
Base.:^(v::AbstractAlgebra.Generic.FreeModuleElem{T},x::MatrixGroupElem{T}) where T <: FieldElem = v.parent(v.v*x.elm)

# evaluation of the form x into the vectors v and u
Base.:*(v::AbstractAlgebra.Generic.FreeModuleElem{T},x::MatrixGroupElem{T},u::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: FieldElem = (v.v*x.elm*transpose(u.v))[1]

map(f::Function, v::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: FieldElem = v.parent(map(f,v.v))
