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

function (a::Generic.RationalFunctionFieldElem)(b::RingElem)
  return divexact(numerator(a)(b), denominator(a)(b))
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

function evaluate(f::AbstractAlgebra.Generic.FracFieldElem{<:MPolyRingElem}, a::Vector{T}) where {T<:RingElem}
  return evaluate(numerator(f), a)//evaluate(denominator(f), a)
end

function evaluate(f::AbstractAlgebra.Generic.FracFieldElem{<:PolyRingElem}, a::RingElem)
  return evaluate(numerator(f), a)//evaluate(denominator(f), a)
end
