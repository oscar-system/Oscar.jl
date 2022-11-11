########################################################################
# Methods for PrincipalOpenSubsets                                     #
########################################################################


########################################################################
# Preimages of PrincipalOpenSubsets                                    #
########################################################################
function preimage(f::AbsSpecMor, U::PrincipalOpenSubset; check::Bool=true) 
  if ambient_scheme(U) != codomain(f) 
    Z = preimage(f, ambient_scheme(U), check=check)
    return PrincipalOpenSubset(Z, OO(Z)(pullback(f)(complement_equation(U)), check=false))
  end
  return PrincipalOpenSubset(domain(f), pullback(f)(complement_equation(U)))
end

########################################################################
# Generic fractions on PrincipalOpenSubsets                            #
########################################################################
@Markdown.doc """
    generic_fraction(a::MPolyLocalizedRingElem, U::PrincipalOpenSubset)

Given a regular function ``a ∈ 𝒪(U)`` on a principal open 
subset ``U ⊂ X`` of an affine scheme ``X``, return a 
fraction ``p/q`` in `Quot(P)` (where ``P`` is the `ambient_coordinate_ring` 
of ``U``) which represents ``a``
in the sense that the maximal extension of its restriction 
to ``U`` returns ``a``.
"""
function generic_fraction(a::MPolyLocalizedRingElem, U::PrincipalOpenSubset)
  X = ambient_scheme(U)
  parent(a) == OO(U) || error("domains are not compatible")
  return lifted_numerator(a)//lifted_denominator(a)
end

function generic_fraction(a::MPolyQuoLocalizedRingElem, U::PrincipalOpenSubset)
  X = ambient_scheme(U)
  parent(a) == OO(U) || error("domains are not compatible")
  return lifted_numerator(a)//lifted_denominator(a)
end

