###############################################################################
# A place to accumulate code that should eventually be moved to AbstractAlgebra.jl
###############################################################################

function Base.copy(f::MPolyRingElem)
  Ox = parent(f)
  g = MPolyBuildCtx(Ox)
  for (c, e) in Base.Iterators.zip(MPolyCoeffs(f), MPolyExponentVectors(f))
    push_term!(g, c, e)
  end
  return finish(g)
end

########################################################################
# Part of PR #4706
function is_equal_as_morphism(f::Map, g::Map)
  f === g && return true
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false
  error("comparison of maps $f and $g not possible")
end
# end of changes in PR #4706
########################################################################


domain_type(::Type{Map{D, C}}) where {D, C} = D
domain_type(f::Map) = domain_type(typeof(f))
codomain_type(::Type{Map{D, C}}) where {D, C} = C
codomain_type(f::Map) = codomain_type(typeof(f))

########################################################################
# Part of PR #5436
@doc raw"""
    matrix(A::Vector{AbstractAlgebra.Generic.FreeModuleElem})

Return the matrix whose rows are the vectors in `A`.
All vectors in `A` must have the same length and the same base ring.

# Examples
```jldoctest
julia> V = vector_space(GF(2), 2); matrix(gens(V))
[1   0]
[0   1]
```
"""
function matrix(A::Vector{AbstractAlgebra.Generic.FreeModuleElem{T}}) where T <: RingElem
   c = length(A[1].v)
   @assert all(x -> length(x.v)==c, A) "Vectors must have the same length"
   X = zero_matrix(base_ring(A[1]), length(A), c)
   for i in 1:length(A), j in 1:c
      X[i,j] = A[i][j]
   end

   return X
end

@doc raw"""
    complement(V::AbstractAlgebra.Generic.FreeModule{T}, W::AbstractAlgebra.Generic.Submodule{T}) where T <: FieldElem

Return a complement for `W` in `V`, i.e. a subspace `U` of `V` such that `V` is
the direct sum of `U` and `W`.

# Examples
```jldoctest
julia> V = vector_space(GF(2), 2);

julia> U = sub(V,[V[1] + V[2]])[1];

julia> C, emb = complement(V,U);

julia> map(emb, gens(C))[1]
(0, 1)
```
"""
function complement(V::AbstractAlgebra.Generic.FreeModule{T}, W::AbstractAlgebra.Generic.Submodule{T}) where T <: FieldElem
   @assert is_submodule(V,W) "The second argument is not a subspace of the first one"
   if vector_space_dim(W)==0 return sub(V,basis(V)) end

   e = W.map

   H = matrix( vcat([e(g) for g in gens(W)], [zero(V) for i in 1:(vector_space_dim(V)-vector_space_dim(W)) ]) )
   A_left = identity_matrix(base_ring(V), vector_space_dim(V))
   A_right = identity_matrix(base_ring(V), vector_space_dim(V))
   for rn in 1:vector_space_dim(W)     # rn = row number
      cn = rn    # column number
      while H[rn,cn]==0 cn+=1 end   # bring on the left the first non-zero entry
      swap_cols!(H,rn,cn)
      swap_rows!(A_right,rn,cn)
      for j in rn+1:vector_space_dim(W)
         add_row!(H,H[j,rn]*H[rn,rn]^-1,rn,j)
         add_column!(A_left,A_left[j,rn]*A_left[rn,rn]^-1,j,rn)
      end
   end
   for j in vector_space_dim(W)+1:vector_space_dim(V)  H[j,j]=1  end
   H = A_left*H*A_right
   _gens = [V([H[i,j] for j in 1:vector_space_dim(V)]) for i in vector_space_dim(W)+1:vector_space_dim(V) ]

   return sub(V,_gens)
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
# end of changes in PR #5436
########################################################################
