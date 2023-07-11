@attributes mutable struct RationalMap{DomainType<:AbsCoveredScheme, 
                                       CodomainType<:AbsCoveredScheme} 
  domain::DomainType
  codomain::CodomainType
  domain_covering::Covering
  codomain_covering::Covering
  domain_chart::AbsSpec
  codomain_chart::AbsSpec
  coord_imgs::Vector{<:FieldElem}

  patch_representatives::IdDict{<:AbsSpec, <:Tuple{<:AbsSpec, <:Vector{<:FieldElem}}}
  realizations::IdDict{<:AbsSpec, <:Vector{<:AbsSpecMor}}

  function RationalMap(
      X::AbsCoveredScheme, Y::AbsCoveredScheme, 
      U::AbsSpec, V::AbsSpec,
      a::Vector{<:FieldElem};
      check::Bool=true,
      domain_covering::Covering=default_covering(X),
      codomain_covering::Covering=default_covering(Y)
    )
    @check is_irreducible(X) "domain must be irreducible"
    @check is_irreducible(Y) "codomain must be irreducible"
    #_find_chart(U, default_covering(X)) !== nothing || error("patch not found in domain")
    #_find_chart(V, default_covering(Y)) !== nothing || error("patch not found in codomain")
    any(x->x===U, patches(default_covering(X))) || error("patch not found in domain")
    any(x->x===V, patches(default_covering(Y))) || error("patch not found in codomain")
    F = parent(first(a))
    R = base_ring(F)
    all(x->parent(x)===F, a) || error("coordinate images must be elements of the same field")
    R === ambient_coordinate_ring(U) || error("images of pullback of the coordinates do not live in the correct ring")
    patch_repr = IdDict{AbsSpec, Tuple{AbsSpec, Vector{FieldElem}}}()
    patch_repr[U] = (V, a)
    realizations = IdDict{AbsSpec, Vector{AbsSpecMor}}()
    return new{typeof(X), typeof(Y)}(X, Y, domain_covering, codomain_covering, 
                                     U, V, a, patch_repr, realizations)
  end
end

domain(Phi::RationalMap) = Phi.domain
codomain(Phi::RationalMap) = Phi.codomain
domain_covering(Phi::RationalMap) = Phi.domain_covering
codomain_covering(Phi::RationalMap) = Phi.codomain_covering
domain_chart(Phi::RationalMap) = Phi.domain_chart
codomain_chart(Phi::RationalMap) = Phi.codomain_chart
coordinate_images(Phi::RationalMap) = Phi.coord_imgs

patch_representatives(Phi::RationalMap) = Phi.patch_representatives
realizations(Phi::RationalMap) = Phi.realizations

function realize_on_patch(Phi::RationalMap, U::AbsSpec)
  if haskey(realizations(Phi), U)
    return realizations(Phi)[U]
  end
  X = domain(Phi)
  Y = codomain(Phi)
  V = codomain_chart(Phi)

  # Try to cover U by PrincipalOpenSubsets W so that the restriction 
  # of Phi to W extends to a regular morphism φ : W → V' for some 
  # `affine_chart` of the codomain of Phi.
  covered_codomain_patches = Vector{AbsSpec}([codomain_chart(Phi)])
  complement_equations = Vector{elem_type(OO(U))}()
  FY = FunctionField(Y)
  FX = FunctionField(X)
  A = [FX(a) for a in coordinate_images(Phi)]
  a = [b[U] for b in A]
  list_for_V = _extend(U, a)
  Psi = [SpecMor(W, ambient_space(V), b, check=false) for (W, b) in list_for_V]
  # Up to now we have maps to the ambient space of V. 
  # But V might be a hypersurface complement in there and we 
  # might need to restrict our domain of definition accordingly. 
  Psi_res = [_restrict_properly(psi, V) for psi in Psi]
  append!(complement_equations, [OO(U)(lifted_numerator(complement_equation(domain(psi)))) for psi in Psi_res])
  while !isone(ideal(OO(U), complement_equations))
    # Find another chart in the codomain which is hopefully easily accessible
    V_next, V_orig = _find_good_neighboring_patch(codomain_covering(Phi), covered_codomain_patches)
    # Get the glueing morphisms for the glueing to some already covered chart
    f, g = glueing_morphisms(glueings(codomain_covering(Phi))[(V_next, V_orig)])
    # Find one morphism which was already realized with this codomomain
    phi = first([psi for psi in Psi_res if codomain(psi) === V_orig])
    # We need to express the pullback of the coordinates of V_next as rational functions 
    # first on V_orig and then pulled back to U
    y0 = gens(OO(V_orig))
    y1 = gens(OO(V_next))
    pb_y1 = pullback(g).(y1)
    rat_lift_y1 = [lifted_numerator(a)//lifted_denominator(a) for a in pb_y1]
    pb_y0 = pullback(phi).(y0)
    rat_lift_y0 = [lifted_numerator(a)//lifted_denominator(a) for a in pb_y0]
    total_rat_lift = [evaluate(a, rat_lift_y0) for a in rat_lift_y1]
    list_for_V_next = _extend(U, total_rat_lift)
    Psi = [SpecMor(W, ambient_space(V_next), b, check=false) for (W, b) in list_for_V_next]
    Psi = [_restrict_properly(psi, V_next) for psi in Psi]
    append!(Psi_res, Psi)
    append!(complement_equations, [OO(U)(lifted_numerator(complement_equation(domain(psi)))) for psi in Psi])
  end
  return Psi_res
end

function realize(Phi::RationalMap)
  realizations = AbsSpecMor[]
  mor_dict = IdDict{AbsSpec, AbsSpecMor}()
  for U in patches(domain_covering(Phi))
    loc_mors = realize_on_patch(Phi, U)
    for phi in loc_mors 
      mor_dict[domain(phi)] = phi
    end
    append!(realizations, loc_mors)
  end
  domain_ref = Covering([domain(phi) for phi in realizations])
  inherit_glueings!(domain_ref, domain_covering(Phi))
  # TODO: Inherit the decomposition_info, too!
  phi_cov = CoveringMorphism(domain_ref, codomain_covering(Phi), mor_dict)
  # Make the refinement known to the domain
  push!(coverings(domain(Phi)), domain_ref)
  return CoveredSchemeMorphism(domain(Phi), codomain(Phi), phi_cov)
end

function _extend(U::AbsSpec, a::Vector{<:FieldElem})
  R = ambient_coordinate_ring(U)
  if iszero(length(a))
    return [(U, elem_type(U)[])]
  end
  F = parent(first(a))
  all(x->parent(x)===F, a) || error("elements must belong to the same field")
  R === base_ring(F) || error("base_rings are incompatible")

  # Determine an ideal for the complement of the maximal domain of definition 
  # for all the a's.
  I_undef = ideal(OO(U), one(OO(U)))
  for f in a
    J = quotient(ideal(OO(U), denominator(f)), ideal(OO(U), numerator(f)))
    I_undef = intersect(I_undef, J)
  end
  #I_undef = ideal(OO(U), small_generating_set(I_undef))

  result = Vector{Tuple{AbsSpec, Vector{RingElem}}}()

  for g in small_generating_set(I_undef)
    Ug = PrincipalOpenSubset(U, g)
    b = [divides(OO(Ug)(numerator(f)), OO(Ug)(denominator(f)))[2] for f in a]
    push!(result, (Ug, b))
  end

  return result
end

function _find_good_neighboring_patch(cov::Covering, covered::Vector{<:AbsSpec})
  U = [x for x in patches(cov) if !any(y->y===x, covered)]
  glue = glueings(cov)
  good_neighbors = [(x, y) for x in U for y in covered if 
                    haskey(glue, (x, y)) && 
                    (glue[(x, y)] isa SimpleGlueing || 
                     (glue[(x, y)] isa LazyGlueing && is_computed(glue[(x, y)]))
                    )
                   ]
  if !isempty(good_neighbors)
    return first(good_neighbors)
  end
  return first(U), first(covered)
end

function _restrict_properly(f::AbsSpecMor, V::AbsSpec{<:Ring, <:MPolyRing})
  return f
end

function _restrict_properly(f::AbsSpecMor, V::AbsSpec{<:Ring, <:MPolyQuoRing})
  return f
end

function _restrict_properly(
    f::AbsSpecMor{<:PrincipalOpenSubset}, V::AbsSpec{<:Ring, <:RT}
  ) where {RT<:MPolyLocRing{<:Ring, <:RingElem, 
                            <:MPolyRing, <:MPolyRingElem, 
                            <:MPolyPowersOfElement}
          }
  h = denominators(inverted_set(OO(V)))
  pbh = pullback(f).(h)
  U = domain(f)
  W = ambient_scheme(f)
  UU = PrincipalOpenSubset(W, push!(pbh, complement_equation(U)))
  return restrict(f, UU, V, check=false)
end

function _restrict_properly(
    f::AbsSpecMor{<:PrincipalOpenSubset}, V::AbsSpec{<:Ring, <:RT}
  ) where {RT<:MPolyQuoLocRing{<:Ring, <:RingElem, 
                            <:MPolyRing, <:MPolyRingElem, 
                            <:MPolyPowersOfElement}
          }
  h = denominators(inverted_set(OO(V)))
  pbh = pullback(f).(h)
  U = domain(f)
  W = ambient_scheme(f)
  UU = PrincipalOpenSubset(W, push!(pbh, complement_equation(U)))
  return restrict(f, UU, V, check=false)
end

