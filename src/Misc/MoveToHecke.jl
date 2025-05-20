###############################################################################
# A place to accumulate code that should eventually be moved to Hecke.jl
###############################################################################

canonical_unit(x::AbsSimpleNumFieldOrderQuoRingElem) = one(parent(x))


########################################################################
# Part of PR #4706
function is_equal_as_morphism(f::MapFromFunc, g::MapFromFunc)
  f === g && return true
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false
  f.f === g.f && return true
  error("comparison of maps $f and $g not possible")
end
# end of changes in PR #4706
########################################################################
