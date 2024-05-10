###############################################################################
# A place to accumulate code that should eventually be moved to AbstractAlgebra.jl
###############################################################################

function symbols(Kt::Generic.RationalFunctionField)
  return Kt.S
end

function number_of_generators(
  F::AbstractAlgebra.Generic.FracField{T}
) where {T<:MPolyRingElem}
  return number_of_generators(base_ring(F))
end

function gen(F::AbstractAlgebra.Generic.FracField{T}) where {T<:PolyRingElem}
  return F(gen(base_ring(F)))
end

function gen(F::AbstractAlgebra.Generic.FracField{T}, i::Int) where {T<:MPolyRingElem}
  return F(gen(base_ring(F), i))
end

function gens(
  F::AbstractAlgebra.Generic.FracField{T}
) where {T<:Union{PolyRingElem,MPolyRingElem}}
  return map(F, gens(base_ring(F)))
end

"""
    is_alternating(B::MatElem)

Return whether the form corresponding to the matrix `B` is alternating,
i.e. `B = -transpose(B)` and `B` has zeros on the diagonal.
Return `false` if `B` is not a square matrix.
"""
function is_alternating(B::MatElem)
  n = nrows(B)
  n == ncols(B) || return false

  for i in 1:n
    B[i, i] == 0 || return false
    for j in (i + 1):n
      B[i, j] == -B[j, i] || return false
    end
  end

  return true
end

@doc raw"""
    upper_triangular_matrix(L)

Return the upper triangular matrix whose entries on and above the diagonal are the elements of `L`.

An exception is thrown whenever the length of `L` is not equal to $n(n+1)/2$,
for some integer $n$.
"""
function upper_triangular_matrix(L)
  T = eltype(L)
  @assert T <: RingElem "L must be a collection of ring elements"
  d = Int(floor((sqrt(1 + 8 * length(L)) - 1) / 2))
  @req length(L) == div(d * (d + 1), 2) "Input vector of invalid length"
  R = parent(L[1])
  x = zero_matrix(R, d, d)
  pos = 1
  for i in 1:d, j in i:d
    x[i, j] = L[pos]
    pos += 1
  end
  return x
end

@doc raw"""
    lower_triangular_matrix(L)

Return the upper triangular matrix whose entries on and below the diagonal are the elements of `L`.

An exception is thrown whenever the length of `L` is not equal to $n(n+1)/2$,
for some integer $n$.
"""
function lower_triangular_matrix(L)
  T = eltype(L)
  @assert T <: RingElem "L must be a collection of ring elements"
  d = Int(floor((sqrt(1 + 8 * length(L)) - 1) / 2))
  @req length(L) == div(d * (d + 1), 2) "Input vector of invalid length"
  R = parent(L[1])
  x = zero_matrix(R, d, d)
  pos = 1
  for i in 1:d, j in 1:i
    x[i, j] = L[pos]
    pos += 1
  end
  return x
end

function Base.copy(f::MPolyRingElem)
  Ox = parent(f)
  g = MPolyBuildCtx(Ox)
  for (c, e) in Base.Iterators.zip(MPolyCoeffs(f), MPolyExponentVectors(f))
    push_term!(g, c, e)
  end
  return finish(g)
end
