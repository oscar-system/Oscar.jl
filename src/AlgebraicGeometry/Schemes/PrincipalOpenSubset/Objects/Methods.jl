########################################################################
# Methods for PrincipalOpenSubsets                                     #
########################################################################


########################################################################
# Preimages of PrincipalOpenSubsets                                    #
########################################################################
function preimage(f::AbsAffineSchemeMor, U::PrincipalOpenSubset; check::Bool=true) 
  if ambient_scheme(U) != codomain(f) 
    Z = preimage(f, ambient_scheme(U), check=check)
    h = lifted_numerator(complement_equation(U))
    fac = factor(h)
    pbh = [OO(Z)(pullback(f)(a)) for (a, _) in fac]
    return PrincipalOpenSubset(Z, prod(pbh; init=one(OO(Z))))
  end
  h = lifted_numerator(complement_equation(U))
  fac = factor(h)
  pbh = [OO(domain(f))(pullback(f)(a)) for (a, _) in fac]
  return PrincipalOpenSubset(domain(f), prod(pbh; init=one(OO(domain(f)))))
end

function preimage(f::AbsAffineSchemeMor{<:AbsAffineScheme, <:PrincipalOpenSubset}, U::PrincipalOpenSubset; check::Bool=true) 
  if ambient_scheme(U) === ambient_scheme(codomain(f))
    h = lifted_numerator(complement_equation(U))
    fac = factor(h)
    pbh = [OO(domain(f))(pullback(f)(a)) for (a, _) in fac]
    return PrincipalOpenSubset(domain(f), prod(pbh; init=one(OO(domain(f)))))
    #return PrincipalOpenSubset(domain(f), pullback(f)(lifted_numerator(complement_equation(U))))
  elseif ambient_scheme(U) === codomain(f)
    h = lifted_numerator(complement_equation(U))
    fac = factor(h)
    pbh = [OO(domain(f))(pullback(f)(a)) for (a, _) in fac]
    return PrincipalOpenSubset(domain(f), prod(pbh; init=one(OO(domain(f)))))
    #return PrincipalOpenSubset(domain(f), pullback(f)(complement_equation(U)))
  end
  # TODO: Make use of the tree structure induced for PrincipalOpenSubset to extend the above pattern.
  Z = preimage(f, ambient_scheme(U), check=check)
  h = lifted_numerator(complement_equation(U))
  fac = factor(h)
  pbh = [OO(Z)(pullback(f)(a)) for (a, _) in fac]
  return PrincipalOpenSubset(Z, prod(pbh; init=one(OO(Z))))
  #return PrincipalOpenSubset(Z, OO(Z)(pullback(f)(lifted_numerator(complement_equation(U)))))
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
  parent(a) == OO(U) || error("domains are not compatible")
  return lifted_numerator(a)//lifted_denominator(a)
end

function generic_fraction(a::MPolyQuoLocRingElem, U::PrincipalOpenSubset)
  parent(a) == OO(U) || error("domains are not compatible")
  return lifted_numerator(a)//lifted_denominator(a)
end

########################################################################
# Base change
########################################################################

function base_change(phi::Any, U::PrincipalOpenSubset;
    ambient_map::AbsAffineSchemeMor=base_change(phi, ambient_scheme(U))[2] # the base change on the ambient scheme
  )
  Y = domain(ambient_map)
  pbf = pullback(ambient_map)
  h = pbf(complement_equation(U))
  UU = PrincipalOpenSubset(Y, h)
  @assert _has_coefficient_map(pullback(ambient_map))
  res_map = restrict(ambient_map, UU, U, check=false)
  #@assert _has_coefficient_map(pullback(res_map))
  return UU, res_map
end

