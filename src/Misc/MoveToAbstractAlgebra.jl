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
    is_zero_entry(B, i, i) || return false
    for j in (i + 1):n
      B[i, j] == -B[j, i] || return false
    end
  end

  return true
end

function Base.copy(f::MPolyRingElem)
  Ox = parent(f)
  g = MPolyBuildCtx(Ox)
  for (c, e) in Base.Iterators.zip(MPolyCoeffs(f), MPolyExponentVectors(f))
    push_term!(g, c, e)
  end
  return finish(g)
end

function (a::Generic.RationalFunctionFieldElem)(b::RingElem)
  return divexact(numerator(a)(b), denominator(a)(b))
end

function evaluate(f::AbstractAlgebra.Generic.FracFieldElem{<:MPolyRingElem}, a::Vector{T}) where {T<:RingElem}
  return evaluate(numerator(f), a)//evaluate(denominator(f), a)
end

function evaluate(f::AbstractAlgebra.Generic.FracFieldElem{<:PolyRingElem}, a::RingElem)
  return evaluate(numerator(f), a)//evaluate(denominator(f), a)
end
