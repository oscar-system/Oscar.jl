@doc raw"""
    RationalMap{DomainType<:AbsCoveredScheme, CodomainType<:AbsCoveredScheme} 

A lazy type for a morphism ``φ : X → Y`` of `AbsCoveredScheme`s which is given 
by a set of rational functions ``a₁,…,aₙ`` in the fraction field of the `base_ring`
of ``𝒪(U)`` for one of the dense open `affine_chart`s ``U`` of ``X``. 
The ``aᵢ`` represent the pullbacks of the coordinates (`gens`) of some 
`affine_chart` ``V`` of the codomain ``Y`` under this map. 
```jldoctest
julia> IP1 = covered_scheme(projective_space(QQ, [:s, :t]))
Scheme
  over rational field
with default covering
  described by patches
    1: spec of multivariate polynomial ring
    2: spec of multivariate polynomial ring
  in the coordinate(s)
    1: [(t//s)]
    2: [(s//t)]

julia> IP2 = covered_scheme(projective_space(QQ, [:x, :y, :z]))
Scheme
  over rational field
with default covering
  described by patches
    1: spec of multivariate polynomial ring
    2: spec of multivariate polynomial ring
    3: spec of multivariate polynomial ring
  in the coordinate(s)
    1: [(y//x), (z//x)]
    2: [(x//y), (z//y)]
    3: [(x//z), (y//z)]

julia> U = first(affine_charts(IP1))
Spectrum
  of multivariate polynomial ring in 1 variable (t//s)
    over rational field

julia> V = first(affine_charts(IP2))
Spectrum
  of multivariate polynomial ring in 2 variables (y//x), (z//x)
    over rational field

julia> t = first(gens(OO(U)))
(t//s)

julia> Phi = oscar.RationalMap(IP1, IP2, U, V, [1//t, 1//t^2]);

julia> realizations = oscar.realize_on_patch(Phi, U);

julia> realizations[3]
Morphism
  from [(t//s)]          spec of localized ring
  to   [(x//z), (y//z)]  spec of multivariate polynomial ring
given by the pullback function
  (x//z) -> (t//s)^2
  (y//z) -> (t//s)

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
  realization_previews::IdDict{<:Tuple{<:AbsSpec, <:AbsSpec}, <:Vector{<:FieldElem}}
  maximal_extensions::IdDict{<:Tuple{<:AbsSpec, <:AbsSpec}, <:Vector{<:AbsSpecMor}}
  cheap_realizations::IdDict{<:Tuple{<:AbsSpec, <:AbsSpec}, <:AbsSpecMor}
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
    realization_previews = IdDict{Tuple{AbsSpec, AbsSpec}, Vector{FieldElem}}()
    maximal_extensions = IdDict{Tuple{AbsSpec, AbsSpec}, Vector{AbsSpecMor}}()
    cheap_realizations = IdDict{Tuple{AbsSpec, AbsSpec}, AbsSpecMor}()
    return new{typeof(X), typeof(Y)}(X, Y, domain_covering, codomain_covering, 
                                     U, V, a, patch_repr, realizations, 
                                     realization_previews, maximal_extensions,
                                     cheap_realizations
                                    )
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
maximal_extensions(Phi::RationalMap) = Phi.maximal_extensions
realization_previews(Phi::RationalMap) = Phi.realization_previews
cheap_realizations(Phi::RationalMap) = Phi.cheap_realizations

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
  #a = [lift(simplify(OO(U)(numerator(b))))//lift(simplify(OO(U)(denominator(b)))) for b in a]
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
    #total_rat_lift = [lift(simplify(OO(U)(numerator(b))))//lift(simplify(OO(U)(denominator(b)))) for b in total_rat_lift]
    list_for_V_next = _extend(U, total_rat_lift)
    Psi = [SpecMor(W, ambient_space(V_next), b, check=false) for (W, b) in list_for_V_next]
    Psi = [_restrict_properly(psi, V_next) for psi in Psi]
    append!(Psi_res, Psi)
    append!(complement_equations, [OO(U)(lifted_numerator(complement_equation(domain(psi)))) for psi in Psi])
    push!(covered_codomain_patches, V_next)
  end
  realizations(Phi)[U] = Psi_res
  return Psi_res
end

function realize_on_open_subset(Phi::RationalMap, U::AbsSpec, V::AbsSpec)
  X = domain(Phi)
  Y = codomain(Phi)
  # Check that the input is admissible
  if !any(x->x===U, patches(X)) 
    UU = _find_chart(U, default_covering(X))
  end
  if !any(x->x===V, patches(Y)) 
    VV = _find_chart(V, default_covering(Y))
  end
  y = function_field(Y).(gens(OO(V)))
  dom_rep = domain_chart(Phi)
  cod_rep = codomain_chart(Phi)
  y_cod = [a[cod_rep] for a in y]::Vector{<:FieldElem}
  x_dom = [evaluate(a, coordinate_images(Phi)) for a in y_cod]::Vector{<:FieldElem}
  x = function_field(X).(x_dom)
  img_gens_frac = [a[U] for a in x]
  dens = [denominator(a) for a in img_gens_frac]
  U_sub = PrincipalOpenSubset(U, OO(U).(dens))
  img_gens = [OO(U_sub)(numerator(a), denominator(a)) for a in img_gens_frac]
  prelim = SpecMor(U_sub, ambient_space(V), img_gens, check=false) # TODO: Set to false
  return _restrict_properly(prelim, V)
end

function realization_preview(Phi::RationalMap, U::AbsSpec, V::AbsSpec)
  if haskey(realization_previews(Phi), (U, V))
    return realization_previews(Phi)[(U, V)]
  end
  X = domain(Phi)
  Y = codomain(Phi)
  # Check that the input is admissible
  if !any(x->x===U, patches(X)) 
    UU = _find_chart(U, default_covering(X))
  end
  if !any(x->x===V, patches(Y)) 
    VV = _find_chart(V, default_covering(Y))
  end
  y = function_field(Y).(gens(OO(V)))
  dom_rep = domain_chart(Phi)
  cod_rep = codomain_chart(Phi)
  y_cod = [a[cod_rep] for a in y]::Vector{<:FieldElem}
  x_dom = [evaluate(a, coordinate_images(Phi)) for a in y_cod]::Vector{<:FieldElem}
  x = function_field(X).(x_dom)
  img_gens_frac = [a[U] for a in x]
  realization_previews(Phi)[(U, V)] = img_gens_frac
  return img_gens_frac
end

function random_realization(Phi::RationalMap, U::AbsSpec, V::AbsSpec)
  img_gens_frac = realization_preview(Phi, U, V)
  U_sub, img_gens = _random_extension(U, img_gens_frac)
  phi = SpecMor(U_sub, ambient_space(V), img_gens, check=true) # Set to false
  return phi
end

function cheap_realization(Phi::RationalMap, U::AbsSpec, V::AbsSpec)
  if haskey(cheap_realizations(Phi), (U, V))
    return cheap_realizations(Phi)[(U, V)]
  end
  img_gens_frac = realization_preview(Phi, U, V)
  # Try to cancel the fractions heuristically
  for (k, f) in enumerate(img_gens_frac)
    a = numerator(f)
    b = denominator(f)
    aa = OO(U)(a)
    new_num = aa
    new_den = one(aa)
    fac = factor(b)
    for (p, e) in fac
      success, q = divides(aa, OO(U)(p))
      while success && e > 0
        aa = q
        e = e - 1
        success, q = divides(aa, OO(U)(p))
      end
      new_den = new_den * p^e
    end
    @assert aa*denominator(f) == unit(fac)*new_den*numerator(f)
    img_gens_frac[k] = inv(unit(fac))*fraction(aa)//fraction(new_den)
  end
  denoms = [denominator(a) for a in img_gens_frac]
  U_sub = PrincipalOpenSubset(U, OO(U).(denoms))
  img_gens = [OO(U_sub)(numerator(a), denominator(a), check=true) for a in img_gens_frac] # Set to false
  phi = SpecMor(U_sub, ambient_space(V), img_gens, check=true) # Set to false
  cheap_realizations(Phi)[(U, V)] = phi
  return phi
end

function realize_maximally_on_open_subset(Phi::RationalMap, U::AbsSpec, V::AbsSpec)
  if haskey(maximal_extensions(Phi), (U, V))
    return maximal_extensions(Phi)[(U, V)]
  end
  img_gens_frac = realization_preview(Phi, U, V)
  extensions = _extend(U, img_gens_frac)
  result = AbsSpecMor[]
  for (U, g) in extensions
    prelim = SpecMor(U, ambient_space(V), g, check=false)
    push!(result, _restrict_properly(prelim, V))
  end
  maximal_extensions(Phi)[(U, V)] = result
  return result
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
    phi_cov = CoveringMorphism(domain_ref, codomain_covering(Phi), mor_dict, check=false)
    # Make the refinement known to the domain
    push!(coverings(domain(Phi)), domain_ref)
    Phi.full_realization = CoveredSchemeMorphism(domain(Phi), codomain(Phi), phi_cov)
  end
  return Phi.full_realization
end

underlying_morphism(Phi::RationalMap) = realize(Phi)

function _random_extension(U::AbsSpec, a::Vector{<:FieldElem})
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
  #I_undef = radical(I_undef)
# @show I_undef
# @show equidimensional_decomposition_radical(saturated_ideal(I_undef))
# I_undef = last(equidimensional_decomposition_radical(I_undef))

  iszero(I_undef) && error("possible domain of definition is empty")
  min_gens = small_generating_set(I_undef)
  kk = coefficient_ring(R)
  g = sum(rand(kk, 1:10)*f for f in min_gens; init=zero(first(min_gens)))
  while iszero(g)
    g = sum(rand(kk, 1:10)*f for f in min_gens; init=zero(first(min_gens)))
  end
  Ug = PrincipalOpenSubset(U, g)
  b = [OO(Ug)(numerator(f), denominator(f)) for f in a]
  return Ug, b
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
# for (k, f) in enumerate(a)
#   @show f
#   isone(denominator(f)) && continue
#   new_den = one(denominator(f))
#   @show new_den
#   @show typeof(new_den)
#   num = OO(U)(numerator(f))
#   @show length(terms(lift(num)))
#   simplify(num)
#   @show length(terms(lift(num)))
#   for (b, e) in factor(denominator(f))
#     @show b, e
#     aa = simplify(OO(U)(b))
#     @show length(terms(b))
#     @show length(terms(lift(aa)))
#     success, q = divides(num, aa)
#     @show success
#     while success && e > 0
#       num = q
#       e = e - 1
#       success, q = divides(num, aa)
#       @show success, q
#     end
#     @show q
#     new_den = new_den*b^e
#   end
#   @assert numerator(a[k])*new_den == denominator(a[k])*num
#   a[k] = fraction(num)//fraction(new_den)
# end
  for f in a
    J = quotient(ideal(OO(U), denominator(f)), ideal(OO(U), numerator(f)))
    I_undef = intersect(I_undef, J)
  end
  #I_undef = ideal(OO(U), small_generating_set(I_undef))
  #I_undef = radical(I_undef)

  result = Vector{Tuple{AbsSpec, Vector{RingElem}}}()

  for g in small_generating_set(I_undef)
    Ug = PrincipalOpenSubset(U, g)
    b = [OO(Ug)(numerator(f), denominator(f)) for f in a]
    #b = [divides(OO(Ug)(numerator(f)), OO(Ug)(denominator(f)))[2] for f in a]
    push!(result, (Ug, b))
  end

  return result
end

equidimensional_decomposition_radical(I::MPolyQuoIdeal) = [ideal(base_ring(I), gens(J)) for J in equidimensional_decomposition_radical(saturated_ideal(I))]
equidimensional_decomposition_radical(I::MPolyLocalizedIdeal) = [ideal(base_ring(I), gens(J)) for J in equidimensional_decomposition_radical(saturated_ideal(I))]
equidimensional_decomposition_radical(I::MPolyQuoLocalizedIdeal) = [ideal(base_ring(I), gens(J)) for J in equidimensional_decomposition_radical(saturated_ideal(I))]

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
  W = ambient_scheme(U)
  UU = PrincipalOpenSubset(W, push!(OO(W).(lifted_numerator.(pbh)), complement_equation(U)))
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
  W = ambient_scheme(U)
  UU = PrincipalOpenSubset(W, push!(OO(W).(lifted_numerator.(pbh)), complement_equation(U)))
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
    real_patches = collect(keys(realizations(Phi)))
    k = findfirst(x->!isone(I(x)), real_patches)
    U = first(affine_charts(X)) # Assign the variable
    if k === nothing
      k = findfirst(x->!isone(I(x)), affine_charts(X))
      k === nothing && error("no affine chart found on which the component was non-trivial")
      U = affine_charts(X)[k]
    else
      U = real_patches[k]
    end
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

function pullback(phi::RationalMap, C::AbsAlgebraicCycle)
  is_isomorphism(phi) || error("method is currently only implemented for isomorphisms")
  X = domain(phi)
  Y = codomain(phi)
  R = coefficient_ring(C)
  comps = IdDict{IdealSheaf, elem_type(R)}()
  for I in components(C)
    @vprint :RationalMap 1 "trying cheap pullback\n"
    pbI = _try_pullback_cheap(phi, I)
    if pbI === nothing
      @vprint :RationalMap 1 "trying randomized pullback\n"
      pbI = _try_randomized_pullback(phi, I)
      if pbI === nothing
        @vprint :RationalMap 1 "trying the full pullback\n"
        pbI = _pullback(phi, I)
      end
    end
    comps[pbI] = C[I]
  end

  return AlgebraicCycle(X, R, comps)
end

function _try_pullback_cheap(phi::RationalMap, I::IdealSheaf)
  X = domain(phi)
  Y = codomain(phi)
  scheme(I) === Y || error("ideal sheaf not defined on the correct scheme")
  # Find a patch in Y on which this component is visible
  all_V = [V for V in affine_charts(Y) if !isone(I(V))]
  function complexity_codomain(V::AbsSpec)
    return sum(total_degree.(lifted_numerator.(gens(I(V)))); init=0)
  end
  sort!(all_V, lt=(x,y)->complexity_codomain(x)<complexity_codomain(y))
  for V in all_V

    # Find a patch in X in which the pullback is visible
    JJ = IdealSheaf(X)
    all_U = copy(affine_charts(X))
    function complexity(U::AbsSpec)
      a = realization_preview(phi, U, V)
      return maximum(vcat([total_degree(numerator(f)) for f in a], [total_degree(denominator(f)) for f in a]))
    end
    sort!(all_U, lt=(x,y)->complexity(x)<complexity(y))

    # First try to get hold of the component via cheap realizations 
    pullbacks = IdDict{AbsSpec, Ideal}()
    for U in all_U
      psi = cheap_realization(phi, U, V)
      U_sub = domain(psi)
      pullbacks[U] = pullback(psi)(saturated_ideal(I(V)))
    end
      #J = pullback(psi)(saturated_ideal(I(V)))
    function new_complexity(U::AbsSpec)
      return sum(total_degree.(lifted_numerator(gens(pullbacks[U]))); init=0)
    end
    sort!(all_U, lt=(x,y)->new_complexity(x)<new_complexity(y))
    for U in all_U
      J = pullbacks[U]
      if !isone(J)
        JJ = IdealSheaf(X, domain(psi), gens(J))
        return JJ
        break
      end
    end
  end
  return nothing
end

function _try_randomized_pullback(phi::RationalMap, I::IdealSheaf)
  X = domain(phi)
  Y = codomain(phi)
  scheme(I) === Y || error("ideal sheaf not defined on the correct scheme")
  # Find a patch in Y on which this component is visible
  all_V = [V for V in affine_charts(Y) if !isone(I(V))]

  min_var = minimum([ngens(OO(V)) for V in all_V])
  all_V = [V for V in all_V if ngens(OO(V)) == min_var]
  deg_bound = minimum([maximum([total_degree(lifted_numerator(g)) for g in gens(I(V))]) for V in all_V])
  all_V = [V for V in all_V if minimum([total_degree(lifted_numerator(g)) for g in gens(I(V))]) == deg_bound]
  V = first(all_V)

  all_U = copy(affine_charts(X))
  function complexity(U::AbsSpec)
    a = realization_preview(phi, U, V)
    return maximum(vcat([total_degree(numerator(f)) for f in a], [total_degree(denominator(f)) for f in a]))
  end
  sort!(all_U, by=complexity)

  for U in all_U
    psi = random_realization(phi, U, V)
    U_sub = domain(psi)

    J = pullback(psi)(saturated_ideal(I(V)))
    if !isone(J)
      JJ = IdealSheaf(X, domain(psi), gens(J))
      return JJ
    end
  end
  return nothing
end

function _pullback(phi::RationalMap, I::IdealSheaf)
  X = domain(phi)
  Y = codomain(phi)
  scheme(I) === Y || error("ideal sheaf not defined on the correct scheme")
  # Find a patch in Y on which this component is visible
  all_V = [V for V in affine_charts(Y) if !isone(I(V))]

  min_var = minimum(ngens(OO(V)) for V in all_V)
  all_V = [V for V in all_V if ngens(OO(V)) == min_var]
  deg_bound = minimum([maximum([total_degree(lifted_numerator(g)) for g in gens(I(V))]) for V in all_V])
  all_V = [V for V in all_V if minimum([total_degree(lifted_numerator(g)) for g in gens(I(V))]) == deg_bound]
  V = first(all_V)

  all_U = copy(affine_charts(X))
  function complexity(U::AbsSpec)
    a = realization_preview(phi, U, V)
    return maximum(vcat([total_degree(numerator(f)) for f in a], [total_degree(denominator(f)) for f in a]))
  end
  sort!(all_U, lt=(x,y)->complexity(x)<complexity(y))

  for U in all_U
    psi_loc = realize_maximally_on_open_subset(phi, U, V)
    # If we are in different components, skip
    length(psi_loc) > 0 || continue
    J = ideal(OO(domain(first(psi_loc))), elem_type(OO(domain(first(psi_loc))))[])
    cod_ideal = ideal(OO(U), elem_type(OO(U))[])
    for (k, psi) in enumerate(psi_loc)
      if dim(cod_ideal) < dim(I)
        break
      end
      J = pullback(psi)(I(V))
      if !isone(J)
        @assert dim(J) == dim(I)
        JJ = IdealSheaf(X, domain(psi), gens(J))
        return JJ
      end
      cod_ideal = cod_ideal + ideal(OO(U), complement_equation(domain(psi)))
    end
  end
  error("ideal sheaf could not be pulled back")
end

function pullback(phi::RationalMap, D::WeilDivisor)
  return WeilDivisor(pullback(phi)(underlying_cycle(D)), check=false) 
end
