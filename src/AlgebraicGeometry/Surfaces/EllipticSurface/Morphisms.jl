########################################################################
# Reduction to positive characteristic
#
# We allow to store a reduction of an elliptic surface `X` to positive
# characteristic. The user needs to know what they're doing here!
#
# The functionality can be made available by specifying a reduction
# map for the `base_ring` (actually a field) of `X` to a field of
# positive characteristic. This can then be stored in `X` via
# `set_good_reduction_map!`. The latter unlocks certain features such
# as computation of intersection numbers in positive characteristic.
########################################################################
function set_good_reduction_map!(X::EllipticSurface, red_map::Map)
  has_attribute(X, :good_reduction_map) && error("reduction map has already been set")
  kk0 = base_ring(X)
  @assert domain(red_map) === kk0
  kkp = codomain(red_map)
  @assert characteristic(kkp) > 0
  set_attribute!(X, :good_reduction_map=>red_map)
end

function get_good_reduction_map(X::EllipticSurface)
  is_zero(characteristic(base_ring(X))) || error("reduction to positive characteristic is only possible from characteristic zero")
  has_attribute(X, :good_reduction_map) || error("no reduction map is available; please set it manually via `set_good_reduction_map!`")
  return get_attribute(X, :good_reduction_map)::Map
end

@attr Tuple{<:AbsCoveredScheme, <:AbsCoveredSchemeMorphism} function raw_good_reduction(X::EllipticSurface)
  red_map = get_good_reduction_map(X)
  X_red, bc_map = base_change(red_map, X)
  set_attribute!(X_red, :is_irreducible=>true)
  set_attribute!(X_red, :is_reduced=>true)
  set_attribute!(X_red, :is_integral=>true)
  set_attribute!(X_red, :is_equidimensional=>true)
  return X_red, bc_map
end

@attr Map function good_reduction_function_fields(X::EllipticSurface)
  red_map = get_good_reduction_map(X)
  E = generic_fiber(X)
  Ft = base_field(E)
  Pt = base_ring(Ft)
  kk = coefficient_ring(Pt)
  kk_red = codomain(red_map)
  Pt_red, _ = polynomial_ring(kk_red, first(symbols(Pt)); cached=false)
  Ft_red = fraction_field(Pt_red)
  Ft_to_Ft_red = MapFromFunc(Ft, Ft_red, fr->Ft_red(map_coefficients(red_map, numerator(fr); parent=Pt_red), map_coefficients(red_map, denominator(fr); parent=Pt_red)))
  return Ft_to_Ft_red
end

@attr EllipticCurve function good_reduction_generic_fiber(X::EllipticSurface)
  red_map = get_good_reduction_map(X)
  E = generic_fiber(X)
  Ft_to_Ft_red = good_reduction_function_fields(X)
  E_red = base_change(Ft_to_Ft_red, E)
  return E_red
end

@attr Vector{<:EllipticCurvePoint} function good_reduction_rational_points(X::EllipticSurface)
  red_map = good_reduction_function_fields(X)
  result = Vector{EllipticCurvePoint}()
  E_red = good_reduction_generic_fiber(X)
  for P in X.MWL # TODO: Do we have a getter for this?
    if is_infinite(P)
      push!(result, infinity(E_red))
      continue
    end
    push!(result, E_red([red_map(P[1]), red_map(P[2])]))
  end
  return result
end

@attr EllipticSurface function good_reduction(X::EllipticSurface)
  red_map = get_good_reduction_map(X)
  E_red = good_reduction_generic_fiber(X)
  mwl_red = good_reduction_rational_points(X)
  X_red = EllipticSurface(E_red, 2, mwl_red; resolution_strategy=X.resolution_strategy)
end

@attr Tuple{<:MorphismFromRationalFunctions, <:MorphismFromRationalFunctions} function identifications_with_raw_good_reduction(X::EllipticSurface)
  X_red = good_reduction(X)
  W_red = weierstrass_chart_on_minimal_model(X_red)
  X_red_raw, red_raw = raw_good_reduction(X)
  red_raw_cov = covering_morphism(red_raw)
  W_red_raw = domain(first(maps_with_given_codomain(red_raw_cov, weierstrass_chart_on_minimal_model(X))))

  R_red = ambient_coordinate_ring(W_red)
  R_red_raw = ambient_coordinate_ring(W_red_raw)
  raw_to_red = morphism_from_rational_functions(X_red_raw, X_red, W_red_raw, W_red, fraction_field(R_red_raw).(gens(R_red_raw)); check=false)
  set_attribute!(raw_to_red, :is_isomorphism=>true)
  red_to_raw = morphism_from_rational_functions(X_red, X_red_raw, W_red, W_red_raw, fraction_field(R_red).(gens(R_red)); check=false)
  set_attribute!(red_to_raw, :is_isomorphism=>true)
  return raw_to_red, red_to_raw
end


@attr IdDict{AbsWeilDivisor, AbsWeilDivisor} function raw_reduction_of_algebraic_lattice(X::EllipticSurface)
  X_red_raw, bc = raw_good_reduction(X)
  basis_ambient, _, _= algebraic_lattice(X)
  return IdDict{AbsWeilDivisor, AbsWeilDivisor}(D=>_reduce_as_prime_divisor(bc, D) for D in basis_ambient)
end

@attr ZZMatrix function good_reduction_algebraic_lattice(X::EllipticSurface)
  div_red = raw_reduction_of_algebraic_lattice(X)
  from, to = identifications_with_raw_good_reduction(X)
  div_red_pf = IdDict{AbsWeilDivisor, AbsWeilDivisor}(D=>pushforward(from, E) for (D, E) in div_red)
  basis, _, _= algebraic_lattice(X)
  X_red = good_reduction(X)
  red_basis, _, _= algebraic_lattice(X_red)
  result = matrix(ZZ, [div_red_pf[D] == E ? one(ZZ) : zero(ZZ) for D in basis, E in red_basis])
  result[1, 1] = one(ZZ) # identify the generic fibers
  return result
end

@attr MorphismFromRationalFunctions function good_reduction(
    f::MorphismFromRationalFunctions{<:EllipticSurface, <:EllipticSurface}
  )
  X = domain(f)
  @assert X === codomain(f) "reduction to positive characteristic is only implemented for automorphisms"
  W = weierstrass_chart_on_minimal_model(X)
  @assert W === domain_chart(f) === codomain_chart(f) "morphism must be defined on the Weierstrass charts"
  img_gens = coordinate_images(f)
  red_map = get_good_reduction_map(X)

  X_red = good_reduction(X)
  W_red = weierstrass_chart_on_minimal_model(X_red)
  R = ambient_coordinate_ring(W_red)
  FR = fraction_field(R)

  psi = fr -> FR(map_coefficients(red_map, numerator(fr); parent=R), map_coefficients(red_map, denominator(fr); parent=R))

  img_gens_red = psi.(img_gens)
  result = morphism_from_rational_functions(X_red, X_red, W_red, W_red, img_gens_red; check=false)
  has_attribute(f, :is_isomorphism) && get_attribute(f, :is_isomorphism)===true && set_attribute!(result, :is_isomorphism=>true)
  return result
end


### Some functions to do custom pullback of divisors along reduction maps.
# We assume that primeness is preserved along the reduction. In particular, the
# user is responsible for this to hold for all cases used!
# They specify the "good reduction" in the end.
function _reduce_as_prime_divisor(bc::AbsCoveredSchemeMorphism, D::AbsWeilDivisor)
  return WeilDivisor(domain(bc), coefficient_ring(D),
                     IdDict{AbsIdealSheaf, elem_type(coefficient_ring(D))}(
                         _reduce_as_prime_divisor(bc, I) => c for (I, c) in coefficient_dict(D)
                       )
                    )
end

function _reduce_as_prime_divisor(bc::AbsCoveredSchemeMorphism, D::EllipticSurfaceSection)
  P = rational_point(D)
  is_infinite(P) && return _reduce_as_prime_divisor(bc, underlying_divisor(D))
  X = codomain(bc)
  @assert parent(P) === generic_fiber(X)
  W = weierstrass_chart_on_minimal_model(X)
  R = ambient_coordinate_ring(W)
  (x, y, t) = gens(R)
  I = ideal(R, R.([evaluate(denominator(P[1]), t)*x-evaluate(numerator(P[1]), t),
                           evaluate(denominator(P[2]), t)*y-evaluate(numerator(P[2]), t)])
           )
  bc_loc = first(maps_with_given_codomain(covering_morphism(bc), W))
  bc_I = pullback(bc_loc)(I)
  @assert is_one(dim(bc_I))
  set_attribute!(bc_I, :is_prime=>true)
  J = PrimeIdealSheafFromChart(domain(bc), domain(bc_loc), bc_I)
  return WeilDivisor(domain(bc), coefficient_ring(D),
                     IdDict{AbsIdealSheaf, elem_type(coefficient_ring(D))}(J => one(coefficient_ring(D)))
                    )
end

function _reduce_as_prime_divisor(bc::AbsCoveredSchemeMorphism, I::AbsIdealSheaf)
  result = pullback(bc, I)
  has_attribute(I, :_self_intersection) && set_attribute!(result, :_self_intersection=>
                                                          (get_attribute(I, :_self_intersection)::Int))
  return result
end

function _reduce_as_prime_divisor(bc::AbsCoveredSchemeMorphism, I::PrimeIdealSheafFromChart)
  U = original_chart(I)
  bc_cov = covering_morphism(bc)
  V = __find_chart(U, codomain(bc_cov))
  IV = I(V)
  bc_loc = first(maps_with_given_codomain(bc_cov, V))
  J = pullback(bc_loc)(IV)
  set_attribute!(J, :is_prime=>true)
  return PrimeIdealSheafFromChart(domain(bc), domain(bc_loc), J)
end

################################################################################################
#
# Given a rational map f:X --> Y compute f_*: NS(X) -> NS(Y)
#
################################################################################################

function _local_pushforward(loc_map::AbsAffineSchemeMor, I::Ideal)
  U_sub = domain(loc_map)
  E, inc_E = sub(U_sub, I) # The subscheme of the divisor
  E_simp = simplify(E) # Eliminate superfluous variables
  id, id_inv = identification_maps(E_simp)

  comp = compose(compose(id, inc_E), loc_map)

  pb = pullback(comp)
  K = kernel(pb)
  return K
end

function _pushforward_lattice_along_isomorphism(step::MorphismFromRationalFunctions{<:EllipticSurface, <:EllipticSurface})
  @assert is_isomorphism(step) "morphism must be an isomorphism"
  X = domain(step)
  Y = codomain(step)
  UX = weierstrass_chart_on_minimal_model(X)
  UY = weierstrass_chart_on_minimal_model(Y)
  @assert codomain_chart(step) === UY
  fracs = coordinate_images(step)

  WY, _ = weierstrass_model(Y)
  UWY = weierstrass_chart(Y)

  to_weierstrass_Y = morphism_from_rational_functions(X, WY, UX, UWY, fracs, check=false)

  fibration_proj_Y = fibration(Y)

  BY = codomain(fibration_proj_Y)
  UBY = codomain(covering_morphism(fibration_proj_Y)[UY])

  composit = morphism_from_rational_functions(X, BY, UX, UBY, [fracs[3]], check=false)

  lat_X = algebraic_lattice(X)[1]
  if !has_attribute(lat_X[1], :is_prime)
    ex, pt, F = irreducible_fiber(X)
    ex || error("no irreducible fiber found; case not implemented")
    lat_X[1] = weil_divisor(F)
    set_attribute!(lat_X[1], :is_prime=>true)
    set_attribute!(first(components(lat_X[1])), :is_prime=>true)
  end

  # We first estimate for every element in the lattic of X whether its image
  # will be a fiber component, or a (multi-)section.
  pre_select = IdDict{AbsWeilDivisor, AbsIdealSheaf}()

  for D in lat_X
    @assert length(components(D)) == 1 "divisors in the algebraic lattice must be prime"
    I = first(components(D))
    @assert has_is_prime(I) && is_prime(I) "ideal sheaf must be known to be prime"
    pre_select[D] = _pushforward_prime_divisor(composit, I)
  end


  # Now we map them one by one using the knowledge gained above
  result = IdDict{AbsWeilDivisor, AbsWeilDivisor}()
  co_ring = coefficient_ring(zero_section(Y))

  n = length(lat_X)
  mwr = rank(mordell_weil_sublattice(X))
  for (i, D) in enumerate(lat_X)
    @vprint :EllipticSurface 2 "$((i, D, pre_select[D]))\n"
    # D is a non-section
    Q = pre_select[D]
    I = first(components(D))
    @vprint :EllipticSurface 2 "$(typeof(I))\n"
    dom_chart = _find_good_representative_chart(I)
    if i > n - mwr # if this is a section
      dom_chart = weierstrass_chart_on_minimal_model(X)
    end

    if dim(Q) == 0
      @vprint :EllipticSurface 3 "image will be a fiber component\n"
      # find the fiber
      if is_one(Q(UBY)) # fiber over infinity
        # collect all components
        comps = AbsWeilDivisor[]
        for (pt, _, F, E, _) in reducible_fibers(Y)
          if is_zero(pt[2]) # if this is in the fiber over the point at ∞ ∈ ℙ¹
            append!(comps, E[2:end])
          end
        end
        @vprint :EllipticSurface 3 "found total of $(length(comps)) possible components\n"

        # collect all charts
        codomain_charts = AbsAffineScheme[]
        if is_empty(comps) # The fiber over infinity
          @vprint :EllipticSurface 3 "the image must be the fiber over infinity"
          codomain_charts = affine_charts(Y) # TODO: How can we restrict the charts then?
        else
          codomain_charts = AbsAffineScheme[V for V in affine_charts(Y) if any(D->!isone(first(components(D))(V)), comps)]
        end
        @vprint :EllipticSurface 3 "found $(length(codomain_charts)) charts where these components are visible"

        if i > n - mwr # If D is a section
          @vprint :EllipticSurface 3 "divisor to be mapped is a section\n"
          pt = X.MWL[i-(n-mwr)]
          res = _pushforward_section(step, pt; divisor=D, codomain_charts)
          result[D] = WeilDivisor(Y, co_ring, IdDict{AbsIdealSheaf, elem_type(co_ring)}(res::AbsIdealSheaf => one(co_ring)); check=false)
        else
          @vprint :EllipticSurface 3 "divisor to be mapped is NOT a section\n"
          loc_map, dom_chart, cod_chart = _prepare_pushforward_prime_divisor(step, I; domain_chart = dom_chart, codomain_charts)

          loc_map === nothing && error("pushforward preparation did not succeed")
          @assert !is_one(I(domain(loc_map)))
          K = _local_pushforward(loc_map, I(domain(loc_map)))
          @assert !is_one(K)

          JJ = ideal(OO(cod_chart), gens(K))
          res = PrimeIdealSheafFromChart(Y, cod_chart, JJ)

          result[D] = WeilDivisor(Y, co_ring, IdDict{AbsIdealSheaf, elem_type(co_ring)}(res::AbsIdealSheaf => one(co_ring)); check=false)
        end
        continue
      end
      @vprint :EllipticSurface 3 "image will not be in the fiber over infinity\n"

      # fiber over some point ≂̸ ∞.
      t = first(gens(OO(UBY)))

      codomain_charts = copy(affine_charts(Y))

      # Restrict the codomain charts if applicable
      for (i, (p, _, F, E, _)) in enumerate(reducible_fibers(Y))
        p[2] == 0 && continue # Fiber over infinity already caught above
        t0 = p[1]//p[2]
        ideal(OO(UBY), t - t0) == Q(UBY) || continue

        # Collect all patches
        codomain_charts = AbsAffineScheme[V for V in affine_charts(Y) if any(I->!isone(I(V)), components(F))]
        break
      end
      @vprint :EllipticSurface 3 "found $(length(codomain_charts)) charts where these components are visible\n"

      if i > n - mwr # If D is a section
        @vprint :EllipticSurface 3 "divisor to be mapped is a section\n"
        pt = X.MWL[i-(n-mwr)]
        res = _pushforward_section(step, pt; divisor=D, codomain_charts)
        result[D] = WeilDivisor(Y, co_ring, IdDict{AbsIdealSheaf, elem_type(co_ring)}(res::AbsIdealSheaf => one(co_ring)); check=false)
      else
        @vprint :EllipticSurface 3 "divisor to be mapped is NOT a section\n"
        loc_map, dom_chart, cod_chart = _prepare_pushforward_prime_divisor(step, I; codomain_charts)
        loc_map === nothing && error("preparation for pushforward did not succeed")

        @assert !is_one(I(domain(loc_map)))
        K = _local_pushforward(loc_map, I(domain(loc_map)))
        @assert !is_one(K)
        JJ = ideal(OO(cod_chart), gens(K))
        res = PrimeIdealSheafFromChart(Y, cod_chart, JJ)

        result[D] = WeilDivisor(Y, co_ring, IdDict{AbsIdealSheaf, elem_type(co_ring)}(res::AbsIdealSheaf => one(co_ring)); check=false)
      end
    else
      # "pushforward will be a section"
      if i > n - mwr # If D is a section
        pt = X.MWL[i-(n-mwr)]
        res = _pushforward_section(step, pt; divisor=D, codomain_charts=[weierstrass_chart_on_minimal_model(Y)])
        if res === nothing
          # The only section not visible in the weierstrass chart is the zero section
          result[D] = zero_section(Y)
          continue
        end

        result[D] = WeilDivisor(Y, co_ring, IdDict{AbsIdealSheaf, elem_type(co_ring)}(res::AbsIdealSheaf => one(co_ring)); check=false)
      else
        loc_map, dom_chart, cod_chart = _prepare_pushforward_prime_divisor(step, I, domain_chart = dom_chart, codomain_charts = [weierstrass_chart_on_minimal_model(Y)])

        if loc_map === nothing
          # The only section not visible in the weierstrass chart is the zero section
          result[D] = zero_section(Y)
          continue
        end

        @assert !is_one(I(domain(loc_map)))
        K = _local_pushforward(loc_map, I(domain(loc_map)))
        @assert !is_one(K)
        JJ = ideal(OO(cod_chart), gens(K))
        res = PrimeIdealSheafFromChart(Y, cod_chart, JJ)

        result[D] = WeilDivisor(Y, co_ring, IdDict{AbsIdealSheaf, elem_type(co_ring)}(res::AbsIdealSheaf => one(co_ring)); check=false)
      end
    end
  end

  res = AbsWeilDivisor[result[D] for D in lat_X]
  for a in res
    set_attribute!(first(components(a)), :_self_intersection, -2)
  end
  # the first one is the class of the fiber; set that one back
  set_attribute!(first(components(first(res))), :_self_intersection, 0)
  return res
end

function _pushforward_section(
    phi::MorphismFromRationalFunctions{<:EllipticSurface, <:EllipticSurface},
    P::EllipticCurvePoint;
    divisor::AbsWeilDivisor=EllipticSurfaceSection(domain(phi), P),
    codomain_charts::Vector{<:AbsAffineScheme} = affine_charts(codomain(phi))
  )
  X = domain(phi)::EllipticSurface
  Y = codomain(phi)::EllipticSurface
  D = divisor
  I = first(components(D))
  iso, inc = morphism_from_section(X, P; divisor=D)
  U = weierstrass_chart_on_minimal_model(X)
  inc_loc = first(maps_with_given_codomain(inc, U))
  U_C = domain(inc_loc)
  phi_loc, _, V = _prepare_pushforward_prime_divisor(phi, I; domain_chart=U, codomain_charts)
  phi_loc === nothing && return nothing # Indicate that the given selection of codomain charts did not lead to a result
  W = codomain(fibration(X)[U])
  iso_loc = _restrict_properly(cheap_realization(iso, W, U_C), U_C)
  inc_dom_phi_loc = inclusion_morphism(domain(phi_loc))
  UU, to_U_C, to_U = fiber_product(inc_loc, inc_dom_phi_loc)
  WW, a, b = fiber_product(iso_loc, to_U_C)
  psi_loc = compose(compose(b, to_U), phi_loc)
  K = kernel(pullback(psi_loc))
  J = ideal(OO(V), gens(K))
  JJ = PrimeIdealSheafFromChart(Y, V, J)
  return JJ
end

@doc raw"""
    pushforward_on_algebraic_lattices(f::MorphismFromRationalFunctions{<:EllipticSurface, <:EllipticSurface}; algorithm=:default)) -> AbstractSpaceMor

Return the pushforward ``f_*: V_1 \to V_2`` where ``V_i`` is the ambient quadratic space of the `algebraic_lattice`.

This assumes that the image ``f_*(V_1)`` is contained in ``V_2``. If this is not the case, you will get
``f_*`` composed with the orthogonal projection to ``V_2``.

# Algorithm
If the attribute `good_reduction_map` has been set via the internal method `Oscar.set_good_reduction_map!`
then the surfaces and the automorphism can be specialized and the computation carried out after specialization.
This is much faster, especially when working over number fields and for complicated maps `f`.

# Input
The keyword argument `algorithm` can be
- `:default` -- use specialization if possible
- `:specialization` -- use specialization and error if this is not possible
- none of the above -- no specialization
"""
function pushforward_on_algebraic_lattices(f::MorphismFromRationalFunctions{<:EllipticSurface, <:EllipticSurface}; algorithm=:default)
  X1 = domain(f)
  X2 = codomain(f)
  can_specialize =  has_attribute(X1, :good_reduction_map) && has_attribute(X2, :good_reduction_map)
  if algorithm == :specialization || (algorithm==:default && can_specialize)
    match1 = good_reduction_algebraic_lattice(X1)
    match2 = good_reduction_algebraic_lattice(X2)
    f_red = good_reduction(f)
    imgs_divs_red = _pushforward_lattice_along_isomorphism(f_red)
    M_red = matrix([basis_representation(codomain(f_red), i) for i in imgs_divs_red])
    M = match1 * M_red * inv(match2)
  else
    imgs_divs = _pushforward_lattice_along_isomorphism(f)
    M = matrix([basis_representation(codomain(f),i) for i in imgs_divs])
  end
  V1 = ambient_space(algebraic_lattice(domain(f))[3])
  V2 = ambient_space(algebraic_lattice(codomain(f))[3])
  # keep the check on since it is simple compared to all the other computations done here
  fstar = hom(V1,V2, M; check=true)
  return fstar
end

# Given an irreducible divisor D on an elliptic surface X, try to extract a point
# on the generic fiber from it. The return value is `nothing` in case this does not succeed.
function point_on_generic_fiber_from_divisor(I::AbsIdealSheaf{<:EllipticSurface}; check::Bool=true)
  X = scheme(I)
  @check dim(I) == 1 "ideal sheaf must be of dimension one"
  return point_on_generic_fiber_from_divisor(WeilDivisor(X, I; check=false); check)
end

function point_on_generic_fiber_from_divisor(D::AbsWeilDivisor{<:EllipticSurface}; check::Bool=true)
  X = ambient_scheme(D)
  E = generic_fiber(X)
  ex, pt, F = irreducible_fiber(X)
  WF = weil_divisor(F)
  # TODO: Also cover this case by considering the class of a reducible fiber?
  !ex && error("no irreducible fiber exists on this algebraic surface")
  @assert length(components(D)) == 1 "divisor must be irreducible"

  I = first(components(D))
  fib = fibration(X)

  # Check a necessary criterion for being a section
  # J = pushforward(fib, I)
  # is_one(dim(J)) || return nothing
  is_zero(intersect(D, WF)) && return nothing
# @check begin
#   J = pushforward(fib, I)
#   is_one(dim(J))
# end "given divisor can not be a section"

  #@check is_one(intersect(D, WF)) "intersection number with irreducible fiber is not one"

  WX = weierstrass_chart_on_minimal_model(X)
  IWX = I(WX)
  is_one(IWX) && return infinity(E) # Point must be the zero section
  R = ambient_coordinate_ring(WX)
  (x, y, t) = gens(R)

  # In case of a multisection do some extra preparation; see below.
  !is_one(intersect(D, WF)) && return point_on_generic_fiber_from_divisor(_prepare_section(D))

  g = gens(groebner_basis(saturated_ideal(IWX), ordering=lex(gens(R))))

  # extract the coefficients for the section
  kkt = base_field(E)

  # First extract the y-coordinate
  i = findfirst(f->(is_zero(degree(f, 1)) && is_one(degree(f, 2))), g)
  i === nothing && return nothing
  #i === nothing && error("no suitable polynomial found to read off point coordinates")
  f = g[i]
  y_coord = one(kkt)
  ev_vals = [zero(kkt), one(kkt), gen(kkt)]
  num = zero(kkt)
  den = zero(kkt)
  for t in terms(f)
    degree(t, 2) == 1 && (den = den - evaluate(t, ev_vals))
    degree(t, 2) == 0 && (num = num + evaluate(t, ev_vals))
  end
  y_coord = num//den

  # Now extract the x-coordinate
  i = findfirst(f->(is_one(degree(f, 1))), g)
  i === nothing && return nothing
  #i === nothing && error("no suitable polynomial found to read off point coordinates")
  f = g[i]
  x_coord = one(kkt)
  ev_vals = [one(kkt), y_coord, gen(kkt)]
  num = zero(kkt)
  den = zero(kkt)
  for t in terms(f)
    degree(t, 1) == 1 && (den = den - evaluate(t, ev_vals))
    degree(t, 1) == 0 && (num = num + evaluate(t, ev_vals))
  end
  x_coord = num//den

  is_zero(evaluate(equation(E), [x_coord, y_coord])) || return nothing
  #@assert is_zero(evaluate(equation(E), [x_coord, y_coord])) "esteemed point does not lie on the curve"
  P = E([x_coord, y_coord])
  return P
end

# Given an isomorphism phi : X -> Y of elliptic surfaces and a full algebraic lattice L on X,
# push forward the divisors D from L to Y and try to extract points on the generic fiber from
# them.
#
# This returns a list consisting of the points on the generic fiber.
function extract_mordell_weil_basis(phi::MorphismFromRationalFunctions{<:EllipticSurface, <:EllipticSurface})
  X = domain(phi)
  Y = codomain(phi)
  is_isomorphism(phi) || error("morphism must be an isomorphism")
  pf_lat = _pushforward_lattice_along_isomorphism(phi)
  points = EllipticCurvePoint[]
  for D in pf_lat
    P = point_on_generic_fiber_from_divisor(D)
    P === nothing && continue
    push!(points, P)
  end
  return points
end

function _prepare_section(D::AbsWeilDivisor{<:EllipticSurface})
  X = ambient_scheme(D)
  WX = weierstrass_chart_on_minimal_model(X)
  R = ambient_coordinate_ring(WX)
  I = first(components(D))
  IWX = I(WX)
  # We have a multisection in this case.
  # To get a section from it, apply arXiv:2103.15101, Algorithm 1.

  # Build up a helper ring
  kkt = base_field(generic_fiber(X))
  f = equation(generic_fiber(X))
  kktXY = parent(f)
  (xx, yy) = gens(kktXY)
  for (c, e) in zip(coefficients(f), exponents(f))
    if e == [0, 2]
      @assert is_one(c) "polynomial is not normalized"
    end
  end
  f = yy^2 - f # prepare the f from the Lemma
  #kktXY, (xx, yy) = polynomial_ring(kkt, [:X, :Y]; cached=false)

  @assert coefficient_ring(R) === coefficient_ring(base_ring(kkt))
  help_map = hom(R, kktXY, [xx, yy, kktXY(gen(kkt))])

  J = ideal(kktXY, help_map.(gens(saturated_ideal(IWX))))

  J_gens = gens(groebner_basis(J, ordering=lex([yy, xx])))
  i = findfirst(f->degree(f, 2) == 0, J_gens)
  i === nothing && error("assertion of Lemma could not be verified")
  g = J_gens[i]
  i = findfirst(f->degree(f, 2) == 1, J_gens)
  i === nothing && error("assertion of Lemma could not be verified")
  h = J_gens[i]
  c = zero(kkt)
  for t in terms(h)
    if degree(t, 2) == 1
      c = c + evaluate(t, [zero(kkt), one(kkt)])
    end
  end
  !isone(c) && (h = inv(c)*h)
  h = yy - h
  @assert J == ideal(kktXY, [g, yy-h])
  ff = equation(kktXY, generic_fiber(X))
  @assert parent(ff) === parent(f)
  @assert ff == yy^2 - f
  while total_degree(g) > 1
    g = divexact(h^2 - f, g)
    p, q = divrem(h, g)
    h = q
  end

  F = fraction_field(R)
  help_map_back = hom(kktXY, F, u->evaluate(u, F(R[3])), F.([R[1], R[2]]))
  new_gens = [help_map_back(g), help_map_back(yy - h)]
  sec_ideal = ideal(OO(WX), numerator.(new_gens))
  @assert dim(sec_ideal) == 1
  @assert is_prime(sec_ideal)

  # overwrite the local variables
  I = PrimeIdealSheafFromChart(X, WX, sec_ideal)
  return weil_divisor(I)
end

#=
# The map is not dominant and can hence not be realized as a MorphismFromRationalFunctions.
# We keep the code for the moment as it will probably help us to reconstruct this map as a
# proper CoveredSchemeMorphism, once this is needed.
=#
function morphism_from_section(
    X::EllipticSurface, P::EllipticCurvePoint;
    divisor::AbsWeilDivisor=EllipticSurfaceSection(X, P)
  )
  U = weierstrass_chart_on_minimal_model(X)
  II = first(components(divisor))

  # For the zero section we can not use the Weierstrass chart
  if is_infinite(P)
    return id_hom(X)
  end
  @assert !is_one(II(U))

  C, inc_C = sub(II)

  UC = domain(first(maps_with_given_codomain(inc_C, U)))

  B = codomain(fibration(X))
  V = codomain(fibration(X)[weierstrass_chart_on_minimal_model(X)])

  kkt = OO(V)::MPolyRing
  @assert ngens(kkt) == 1
  t = first(gens(kkt))
  img_gens = [evaluate(P.coordx, t), evaluate(P.coordy, t), t]

  Fkkt = fraction_field(kkt)
  img_gens2 = Fkkt.(img_gens)
  # TODO: Cache?
  iso = morphism_from_rational_functions(B, C, V, UC, img_gens2, check=false)
  return iso, inc_C
end

########################################################################
# Translations by sections                                             #
########################################################################

@doc raw"""
    translation_morphism(X::EllipticSurface, P::EllipticCurvePoint) -> MorphismFromRationalFunctions

Return the automorphism of ``X`` defined by fiberwise translation by the section ``P``.
"""
function translation_morphism(X::EllipticSurface, P::EllipticCurvePoint)
  E = generic_fiber(X)
  @assert parent(P) === E "point does not lay on the underlying elliptic curve"
  U = weierstrass_chart_on_minimal_model(X)
  is_zero(P) && return id_hom(X)

  # We construct the translation by P as a morphism of rational functions
  kT = base_field(E)
  T = first(gens(kT))

  R = ambient_coordinate_ring(U)
  x, y, t = gens(R)

  a1, a2, a3, a4, a6 = [evaluate(a, t) for a in a_invariants(E)]

  p_x = evaluate(P[1], t)
  p_y = evaluate(P[2], t)

  # Formulas adapted from Hecke/src/EllCrv/EllCrv.jl
  m = (p_y - y)//(p_x - x)
  pb_x = - x - p_x - a2 + a1*m + m^2
  pb_y = - y - m*(pb_x - x) - a1*pb_x - a3

  F = fraction_field(R)

  result = morphism_from_rational_functions(X, X, U, U, F.([pb_x, pb_y, t]), check=true)
  set_attribute!(result, :is_isomorphism=>true)
  return result
end

########################################################################
# Möbius transformations                                               #
########################################################################

# Find a moebius transformation which sends a given set of three points in ℙ¹ to another set
# of three points.
function find_moebius_transformation(
    orig_pts::Vector{<:Vector{<:FieldElem}},
    new_pts::Vector{<:Vector{<:FieldElem}}
  )
  kk = parent(first(orig_pts))
  a = [a[1] for a in orig_pts]
  b = [b[1] for b in new_pts]
  @assert all(a->isone(a[2]), orig_pts) "not implemented for non-normalized or infinite points"
  @assert all(a->isone(a[2]), new_pts) "not implemented for non-normalized or infinite points"
  return find_moebius_transformation(a, b)
end

function find_moebius_transformation(
    orig_pts::Vector{<:FieldElem},
    new_pts::Vector{<:FieldElem}
  )
  length(orig_pts) == 3 || error("exactly three points are needed")
  @assert length(orig_pts) == length(new_pts) "number of points must coincide"
  kk = parent(first(orig_pts))
  a = orig_pts
  b = new_pts

  # Set up the matrix mapping the first three points to 0, 1, ∞
  A = kk[(a[2] - a[3]) (-a[1]*(a[2] - a[3])); (a[2] - a[1]) (-a[3]*(a[2] - a[1]))]

  # Set up the matrix mapping the second three points to 0, 1, ∞
  B = kk[(b[2] - b[3]) (-b[1]*(b[2] - b[3])); (b[2] - b[1]) (-b[3]*(b[2] - b[1]))]

  C = inv(B)*A
  return x->(C[1,1]*x + C[1, 2], C[2,1]*x + C[2,2])
end

# Given two abstractly isomorphic elliptic surfaces X and Y over ℙ¹,
# find all moebius transformation of the base which preserve the critical
# values of the projections, try to lift them to morphisms X -> Y and
# return the list of such morphisms for which the lift was successful.
@doc raw"""
    isomorphisms(X::EllipticSurface, Y::EllipticSurface) -> Vector{MorphismFromRationalFunctions}

Given two elliptic surfaces `` X \to \mathbb{P}^1`` and `` Y \to \mathbb{P}^1`` return all
isomorphisms ``X \to Y`` such that there exists Möbius transformation ``\mathbb{P}^1 \to \mathbb{P}^1``
fitting in the following commutative diagram.
```math
\begin{array}{ccc}
  X & \to & Y \\
  \downarrow & & \downarrow\\
  \mathbb{P}^1 & \to &  \mathbb{P}^1
\end{array}
```
"""
isomorphisms(X::EllipticSurface, Y::EllipticSurface) = admissible_moebius_transformations(X, Y)

isomorphisms_on_weierstrass_chart(X::EllipticSurface, Y::EllipticSurface) = admissible_moebius_transformations_on_weierstrass_chart(X, Y)

function admissible_moebius_transformations(
    X::EllipticSurface,
    Y::EllipticSurface
  )
  result = MorphismFromRationalFunctions[]
  for img_gens in  _admissible_moebius_transformations(X, Y; on_weierstrass_model=false)
    push!(result, _moebius_to_morphism_from_rational_functions(X, Y, img_gens))
  end
  return result
end

function admissible_moebius_transformations_on_weierstrass_chart(
    X::EllipticSurface,
    Y::EllipticSurface
  )
  result = MapFromFunc[]
  for img_gens in  _admissible_moebius_transformations(X, Y; on_weierstrass_model=true)
    push!(result, _moebius_to_pullback_on_weierstrass_chart(X, Y, img_gens))
  end
  return result
end

function _admissible_moebius_transformations(
    X::EllipticSurface,
    Y::EllipticSurface; on_weierstrass_model=true
  )
  EX = generic_fiber(X)
  EY = generic_fiber(Y)

#  kkt = base_field(EX)
#  @assert kkt === base_field(EY) "base fields of the generic fibers must coincide"
  kk = base_ring(X)
  @assert kk === base_ring(Y) "elliptic surfaces must be defined over the same field"

  dX = numerator(discriminant(EX))::PolyRingElem
  dY = numerator(discriminant(EY))::PolyRingElem

  vX = roots(dX)
  @assert all(is_one(degree(a)) for (a, k) in factor(dX))  "not all critical values are rational over the given ground field"

  vY = roots(dY)
  @assert all(is_one(degree(a)) for (a, k) in factor(dY))  "not all critical values are rational over the given ground field"

#   for (c, _) in reducible_fibers(X)
#     @assert !is_zero(c[2]) "the case of reducible fibers over the point at infinity is not implemented"
#   end
#   for (c, _) in reducible_fibers(Y)
#     @assert !is_zero(c[2]) "the case of reducible fibers over the point at infinity is not implemented"
#   end

  # Use the first three elements of vX and map them to three elements of vY.
  # Then check whether the resulting transformation preserves everything.

  candidates = []

  @assert length(vX) >= 3 "at least three reducible fibers are needed"
  length(vX) == length(vY) || return candidates # No moebius transformation is possible in this case
  kkt = base_field(EX)
  t = gen(kkt)
  p1 = vX[1:3]
  for i in vY
    for j in vY
      i == j && continue
      for k in vY
        (i == k || j == k) && continue
        p2 = [i, j, k]
        mt = find_moebius_transformation(p1, p2)
        any(is_zero(mt(x)[2]) for x in vX) && continue # reducible fibers over ∞ are not implemented at the moment.
        any(!(mt(x)[1]//mt(x)[2] in vY) for x in vX) && continue # the transformation does not preserve all admissible fibers in this case
        p, q = mt(t)
        img_t = (p//q)::typeof(t)
        EYbc = base_change(f->evaluate(f, img_t), EY)
        is_isomorphic(EYbc, EX) || continue
        iso_ell = isomorphism(EX, EYbc)
        push!(candidates, _to_weierstrass_morphism(X, Y, mt, iso_ell; on_weierstrass_model))
      end
    end
  end
  return candidates
end

function _to_weierstrass_morphism(X, Y, mt, iso_ell; on_weierstrass_model)
  EX = generic_fiber(X)
  EY = generic_fiber(Y)
  # Set up some variables
  kkt = base_field(EX)
  t = gen(kkt)
  if on_weierstrass_model
    WX = weierstrass_chart(X)
    WY = weierstrass_chart(Y)
  else
    WX = weierstrass_chart_on_minimal_model(X)
    WY = weierstrass_chart_on_minimal_model(Y)
  end
  RX = ambient_coordinate_ring(WX)
  FRX = fraction_field(RX)
  RY = ambient_coordinate_ring(WY)
  FRY = fraction_field(RY)

  # Construct the isomorphism of elliptic surfaces explicitly

  a, b, _ = rational_maps(iso_ell)
  kkTxy = parent(a)
  to_FRX = hom(kkTxy, FRX, x->evaluate(x, FRX(RX[3])), FRX.([RX[1], RX[2]]))
  A = to_FRX(a)
  B = to_FRX(b)
  P, Q = mt(FRX(RX[3]))
  img_T = (P//Q)::elem_type(FRX)
  img_gens = [A, B, img_T]
  return img_gens
end

function _moebius_to_pullback_on_weierstrass_chart(X, Y, img_gens)
  WY = weierstrass_chart(Y)
  WX = weierstrass_chart(X)
  RX = ambient_coordinate_ring(WX)
  FRX = fraction_field(RX)
  RY = ambient_coordinate_ring(WY)
  FRY = fraction_field(RY)

  return Hecke.extend_domain_to_fraction_field(hom(RY, FRX, img_gens))
end

function _moebius_to_morphism_from_rational_functions(X, Y, img_gens)
  WY = weierstrass_chart_on_minimal_model(Y)
  WX = weierstrass_chart_on_minimal_model(X)

  loc_res = morphism_from_rational_functions(X, Y, WX, WY, img_gens; check=true)
  set_attribute!(loc_res, :is_isomorphism=>true)
  return loc_res
end

# An internal helper routine to verify that a given isomorphism of elliptic surfaces
# does indeed give an isomorphism on their generic fibers.
function check_isomorphism_on_generic_fibers(phi::MorphismFromRationalFunctions{<:EllipticSurface, <:EllipticSurface})
  X = domain(phi)
  Y = codomain(phi)
  @assert domain_chart(phi) === weierstrass_chart_on_minimal_model(X)
  @assert codomain_chart(phi) === weierstrass_chart_on_minimal_model(Y)
  EX = generic_fiber(X)
  EY = generic_fiber(Y)
  a, b, c = coordinate_images(phi)

  hX = equation(EX)
  RX = parent(hX)
  FX = fraction_field(RX)
  kktX = coefficient_ring(RX)

  hY = equation(EY)
  RY = parent(hY)
  FY = fraction_field(RY)
  kktY = coefficient_ring(RY)

  A = evaluate(a, [RX[1], RX[2], RX(gen(kktX))])
  B = evaluate(b, [RX[1], RX[2], RX(gen(kktX))])
  C = evaluate(c, [RX[1], RX[2], RX(gen(kktX))])

  help_map = hom(RY, FX, t->evaluate(t, C), [A, B])

  hh = help_map(hY)

  return divides(hX, numerator(hh))[1]
end

@doc raw"""
    isomorphism_from_generic_fibers(X::EllipticSurface, Y::EllipticSurface, f::Hecke.EllCrvIso) -> MorphismFromRationalFunctions

Given an isomorphism ``f`` between the generic fibers of ``X`` and ``Y``, return the corresponding
isomorphism of ``X`` and ``Y`` over ``\mathbb{P}^1``.
"""
function isomorphism_from_generic_fibers(
    X::EllipticSurface, Y::EllipticSurface, f::Hecke.EllCrvIso
  )
  EX = generic_fiber(X)
  EY = generic_fiber(Y)
  iso_ell = f
  @req domain(f) == EX "must be an isomorphism of the generic fibers"
  @req codomain(f) == EY "must be an isomorphism of the generic fibers"
  a, b, _ = rational_maps(iso_ell)
  kt = base_field(EX)
  t = gen(kt)

  # Make sure we got something reasonable
  h2 = equation(EY)
  pb_h2 = evaluate(h2, [a, b])
  @assert divides(pb_h2, equation(parent(pb_h2), EX))[1]

  WX = weierstrass_chart_on_minimal_model(X)
  RX = ambient_coordinate_ring(WX)
  FRX = fraction_field(RX)
  WY = weierstrass_chart_on_minimal_model(Y)
  RY = ambient_coordinate_ring(WY)
  FRY = fraction_field(RY)

  kkTxy = parent(a)
  to_FRX = hom(kkTxy, FRX, x->evaluate(x, FRX(RX[3])), FRX.([RX[1], RX[2]]))
  A = to_FRX(a)
  B = to_FRX(b)
  img_gens = [A, B, FRX(RX[3])]
  m = morphism_from_rational_functions(X, Y, WX, WY, FRX.(img_gens); check=false)
  set_attribute!(m, :is_isomorphism=>true)
  return m
end

@doc raw"""
    isomorphism_from_generic_fibers(X::EllipticSurface, Y::EllipticSurface) -> MorphismFromRationalFunctions

Given given two elliptic surfaces ``X`` and ``Y`` whose generic fibers are isomorphic,
return the corresponding isomorphism of ``X`` and ``Y`` over ``\mathbb{P}^1``.
"""
function isomorphism_from_generic_fibers(
    X::EllipticSurface, Y::EllipticSurface
  )
  EX = generic_fiber(X)
  EY = generic_fiber(Y)
  is_isomorphic(EX, EY) || error("generic fibers are not isomorphic")
  iso_ell = isomorphism(EX, EY)
  return isomorphism_from_generic_fibers(X, Y, iso_ell)
end
