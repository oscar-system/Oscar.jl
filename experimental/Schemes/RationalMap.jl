@doc raw"""
    RationalMap{DomainType<:AbsCoveredScheme, CodomainType<:AbsCoveredScheme} 

A lazy type for a morphism ``φ : X → Y`` of `AbsCoveredScheme`s which is given 
by a set of rational functions ``a₁,…,aₙ`` in the fraction field of the `base_ring`
of ``𝒪(U)`` for one of the dense open `affine_chart`s ``U`` of ``X``. 
The ``aᵢ`` represent the pullbacks of the coordinates (`gens`) of some 
`affine_chart` ``V`` of the codomain ``Y`` under this map. 
```jldoctest
julia> IP1 = covered_scheme(projective_space(QQ, [:s, :t]))
covered scheme with 2 affine patches in its default covering

julia> IP2 = covered_scheme(projective_space(QQ, [:x, :y, :z]))
covered scheme with 3 affine patches in its default covering

julia> U = first(affine_charts(IP1))
Spec of Multivariate polynomial ring in 1 variable over QQ

julia> V = first(affine_charts(IP2))
Spec of Multivariate polynomial ring in 2 variables over QQ

julia> t = first(gens(OO(U)))
(t//s)

julia> Phi = oscar.RationalMap(IP1, IP2, U, V, [1//t, 1//t^2]);

julia> realizations = oscar.realize_on_patch(Phi, U);

julia> realizations[3]
morphism from

	Spec of Localization of multivariate polynomial ring in 1 variable over QQ at products of 0 elements

to

	Spec of Multivariate polynomial ring in 2 variables over QQ

with coordinates

	(t//s)^2, (t//s)

```
"""
@attributes mutable struct RationalMap{DomainType<:AbsCoveredScheme, 
                                       CodomainType<:AbsCoveredScheme
                                      } <: AbsCoveredSchemeMorphism{DomainType, CodomainType, 
                                                                    RationalMap, Nothing}
  domain::DomainType
  codomain::CodomainType
  domain_covering::Covering
  codomain_covering::Covering
  domain_chart::AbsSpec
  codomain_chart::AbsSpec
  coord_imgs::Vector{<:FieldElem}

  patch_representatives::IdDict{<:AbsSpec, <:Tuple{<:AbsSpec, <:Vector{<:FieldElem}}}
  realizations::IdDict{<:AbsSpec, <:Vector{<:AbsSpecMor}}
  full_realization::CoveredSchemeMorphism

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
  covered_codomain_patches = Vector{AbsSpec}([V])
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
  @assert all(phi->codomain(phi) === V, Psi_res)
  append!(complement_equations, [OO(U)(lifted_numerator(complement_equation(domain(psi)))) for psi in Psi_res])
  while !isone(ideal(OO(U), complement_equations))
    # Find another chart in the codomain which is hopefully easily accessible
    V_next, V_orig = _find_good_neighboring_patch(codomain_covering(Phi), covered_codomain_patches)
    # Get the glueing morphisms for the glueing to some already covered chart
    f, g = glueing_morphisms(glueings(codomain_covering(Phi))[(V_next, V_orig)])
    # Find one morphism which was already realized with this codomomain
    phi = first([psi for psi in Psi_res if codomain(psi) === V_orig])
    # We need to express the pullback of the coordinates of V_next as rational functions, 
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
    push!(covered_codomain_patches, V_next)
  end
  return Psi_res
end

function realize(Phi::RationalMap)
  if !isdefined(Phi, :full_realization)
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
    Phi.full_realization = CoveredSchemeMorphism(domain(Phi), codomain(Phi), phi_cov)
  end
  return Phi.full_realization
end

underlying_morphism(Phi::RationalMap) = realize(Phi)

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
  isempty(U) && error("no new neighbor could be found")
  return first(U), first(covered)
end

function _restrict_properly(f::AbsSpecMor, V::AbsSpec{<:Ring, <:MPolyRing})
  return restrict(f, domain(f), V, check=false)
end

function _restrict_properly(f::AbsSpecMor, V::AbsSpec{<:Ring, <:MPolyQuoRing})
  return restrict(f, domain(f), V, check=false)
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

function pushforward(Phi::RationalMap, D::AbsAlgebraicCycle)
  is_isomorphism(Phi) || error("method not implemented unless for the case of an isomorphism")
  #is_proper(Phi) || error("morphism must be proper")
  all(x->isprime(x), components(D)) || error("divisor must be given in terms of irreducible components")
  X = domain(Phi)
  Y = codomain(Phi)
  pushed_comps = IdDict{IdealSheaf, elem_type(coefficient_ring(D))}()
  for I in components(D)
    # Find some chart in which I is non-trivial
    k = findfirst(x->!isone(I(x)), affine_charts(X))
    k === nothing && error("no affine chart found on which the component was non-trivial")
    U = affine_charts(X)[k]
    loc_phi = realize_on_patch(Phi, U)
    k = findfirst(x->!isone(I(domain(x))), loc_phi)
    k === nothing && error("no patch found on which the component was non-trivial")
    phi = loc_phi[k]
    U = domain(phi)
    V = codomain(phi)
    pb = pullback(phi)
    Q, pr = quo(OO(U), I(U))
    J = kernel(hom(OO(V), Q, compose(pb, pr).(gens(OO(V))), check=false))
    # If this map is contracting the component, skip
    dim(I(U)) == dim(J) || continue
    JJ = IdealSheaf(Y, V, gens(J))
    # TODO: There is a further multiplicity!
    pushed_comps[JJ] = D[I]
  end
  return AlgebraicCycle(Y, coefficient_ring(D), pushed_comps)
end

function pushforward(Phi::RationalMap, D::WeilDivisor)
  return WeilDivisor(pushforward(Phi, underlying_cycle(D)))
end

@attr function is_proper(phi::AbsCoveredSchemeMorphism)
  error("no method implemented to check properness")
end

@attr function is_isomorphism(phi::AbsCoveredSchemeMorphism)
  error("no method implemented to check for being an isomorphism")
end

