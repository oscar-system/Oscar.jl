########################################################################
# Methods for PrincipalOpenSubsets                                     #
########################################################################


########################################################################
# Preimages of PrincipalOpenSubsets                                    #
########################################################################
function preimage(f::AbsSpecMor, U::PrincipalOpenSubset; check::Bool=true) 
  if ambient_scheme(U) != codomain(f) 
    Z = preimage(f, ambient_scheme(U), check=check)
    return PrincipalOpenSubset(Z, OO(Z)(pullback(f)(lifted_numerator(complement_equation(U))), check=false))
  end
  return PrincipalOpenSubset(domain(f), pullback(f)(complement_equation(U)))
end

########################################################################
# Generic fractions on PrincipalOpenSubsets                            #
########################################################################
@doc raw"""
    generic_fraction(a::MPolyLocRingElem, U::PrincipalOpenSubset)

Given a regular function ``a ∈ 𝒪(U)`` on a principal open 
subset ``U ⊂ X`` of an affine scheme ``X``, return a 
fraction ``p/q`` in `Quot(P)` (where ``P`` is the `ambient_coordinate_ring` 
of ``U``) which represents ``a``
in the sense that the maximal extension of its restriction 
to ``U`` returns ``a``.
"""
function generic_fraction(a::MPolyLocRingElem, U::PrincipalOpenSubset)
  X = ambient_scheme(U)
  parent(a) == OO(U) || error("domains are not compatible")
  return lifted_numerator(a)//lifted_denominator(a)
end

function generic_fraction(a::MPolyQuoLocRingElem, U::PrincipalOpenSubset)
  X = ambient_scheme(U)
  parent(a) == OO(U) || error("domains are not compatible")
  return lifted_numerator(a)//lifted_denominator(a)
end

########################################################################
# Base change
########################################################################

function base_change(phi::Any, U::PrincipalOpenSubset;
    ambient_map::AbsSpecMor=base_change(phi, ambient_scheme(U))[2] # the base change on the ambient scheme
  )
  Y = domain(ambient_map)
  pbf = pullback(ambient_map)
  h = pbf(complement_equation(U))
  UU = PrincipalOpenSubset(Y, h)
  return UU, restrict(ambient_map, UU, U, check=true) # TODO: Set to false after testing
end

