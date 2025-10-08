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
