###############################################################################
# A place to accumulate code that should eventually be moved to Nemo.jl
###############################################################################

function minpoly(a::fpFieldElem)
  kx, x = polynomial_ring(parent(a), cached = false)
  return x-a
end

function (F::QQField)(x::AbstractAlgebra.Generic.FracFieldElem{T}) where T <: Union{QQPolyRingElem, QQMPolyRingElem}
  num = numerator(x)
  cst_num = constant_coefficient(num)
  denom = denominator(x)
  cst_denom = constant_coefficient(denom)
  if (num != cst_num || denom != cst_denom) throw(InexactError(:QQFieldElem, QQFieldElem, x)) end
  F(cst_num) // F(cst_denom)
end

function (F::QQField)(x::AbstractAlgebra.Generic.RationalFunctionFieldElem{QQFieldElem})
  return F(x.d)
end

function (F::fpField)(x::AbstractAlgebra.Generic.FracFieldElem{T}) where T <: Union{fpPolyRingElem, fpMPolyRingElem}
  num = numerator(x)
  cst_num = constant_coefficient(num)
  denom = denominator(x)
  cst_denom = constant_coefficient(denom)
  if (num != cst_num || denom != cst_denom) throw(InexactError(:fpFieldElem, fpFieldElem, x)) end
  F(cst_num) // F(cst_denom)
end

function (F::FqField)(x::AbstractAlgebra.Generic.FracFieldElem{T}) where T <: Union{FqPolyRingElem, FqMPolyRingElem}
  num = numerator(x)
  cst_num = constant_coefficient(num)
  denom = denominator(x)
  cst_denom = constant_coefficient(denom)
  if (num != cst_num || denom != cst_denom) throw(InexactError(:FqFieldElem, FqFieldElem, x)) end
  F(cst_num) // F(cst_denom)
end

function (F::fpField)(x::AbstractAlgebra.Generic.RationalFunctionFieldElem{fpFieldElem})
  return F(x.d)
end

function (F::FqField)(x::AbstractAlgebra.Generic.RationalFunctionFieldElem{FqFieldElem})
  return F(x.d)
end

^(a::MatElem, b::ZZRingElem) = Nemo._generic_power(a, b)

########################################################################
# Part of PR #5436
@doc raw"""
    conjugate_transpose(x::MatElem{T}) where T <: FinFieldElem

If the base ring of `x` is `GF(q^2)`, return the matrix `transpose( map ( y -> y^q, x) )`.
An exception is thrown if the base ring does not have even degree.

# Examples
```jldoctest
julia> F, z = finite_field(2, 2); M = matrix(F, [0 z; 0 0])
[0   o]
[0   0]

julia> conjugate_transpose(M)
[    0   0]
[o + 1   0]
```
"""
function conjugate_transpose(x::MatElem{T}) where T <: FinFieldElem
  e = degree(base_ring(x))
  @req iseven(e) "The base ring must have even degree"
  e = div(e, 2)

  y = similar(x, ncols(x), nrows(x))
  for i in 1:ncols(x), j in 1:nrows(x)
    # This code could be *much* faster, by precomputing the Frobenius map
    # once; see also FrobeniusCtx in Hecke (but that does not yet support all
    # finite field types at the time this comment was written).
    # If you need this function to be faster, talk to Claus or Max.
    y[i,j] = frobenius(x[j,i],e)
  end
  return y
end

@doc raw"""
    is_hermitian(B::MatElem{T}) where T <: FinFieldElem

Return whether the matrix `B` is hermitian, i.e. `B = conjugate_transpose(B)`.
Return `false` if `B` is not a square matrix, or the field has not even degree.

# Examples
```jldoctest
julia> F, z = finite_field(2, 2); x = matrix(F, [1 z; 0 1]);

julia> is_hermitian(x)
false

julia> is_hermitian(x + conjugate_transpose(x))
true
```
"""
function is_hermitian(B::MatElem{T}) where T <: FinFieldElem
   n = nrows(B)
   n == ncols(B) || return false
   e = degree(base_ring(B))
   iseven(e) ? e = div(e,2) : return false

   for i in 1:n, j in i:n
      B[i,j] == frobenius(B[j,i],e) || return false
   end

   return true
end
# end of changes in PR #5436
########################################################################
