# Type declaration has been moved to Types.jl

domain(Phi::MorphismFromRationalFunctions) = Phi.domain
codomain(Phi::MorphismFromRationalFunctions) = Phi.codomain
domain_covering(Phi::MorphismFromRationalFunctions) = Phi.domain_covering
codomain_covering(Phi::MorphismFromRationalFunctions) = Phi.codomain_covering
domain_chart(Phi::MorphismFromRationalFunctions) = Phi.domain_chart
codomain_chart(Phi::MorphismFromRationalFunctions) = Phi.codomain_chart
coordinate_images(Phi::MorphismFromRationalFunctions) = Phi.coord_imgs

# user facing constructor
@doc raw"""
    morphism_from_rational_functions(
          X::AbsCoveredScheme, Y::AbsCoveredScheme, 
          U::AbsAffineScheme, V::AbsAffineScheme,
          a::Vector{<:FieldElem};
          check::Bool=true,
          domain_covering::Covering=default_covering(X),
          codomain_covering::Covering=default_covering(Y)
        )

Given two `AbsCoveredScheme`s `X` and `Y` this constructs a morphism 
`f : X → Y` from a list of rational functions `a`. The latter are 
interpreted as representatives of the pullback along `f` of the 
`coordinates` of an `affine_chart` `V` of the codomain `Y` in the 
chart `U` in the domain `X`. 

Note that, since there is no type supporting fraction fields of 
quotient rings at the moment, the entries of `a` need to be 
fractions of polynomials in the `ambient_coordinate_ring` of `U`.

```jldoctest
julia> IP1 = covered_scheme(projective_space(QQ, [:s, :t]))
Scheme
  over rational field
with default covering
  described by patches
    1: affine 1-space
    2: affine 1-space
  in the coordinate(s)
    1: [(t//s)]
    2: [(s//t)]

julia> IP2 = projective_space(QQ, [:x, :y, :z]);

julia> S = homogeneous_coordinate_ring(IP2);

julia> x, y, z = gens(S);

julia> IPC, inc_IPC = sub(IP2, ideal(S, [x^2 - y*z]));

julia> C = covered_scheme(IPC);

julia> U = first(affine_charts(IP1))
Spectrum
  of multivariate polynomial ring in 1 variable (t//s)
    over rational field

julia> V = first(affine_charts(C))
Spectrum
  of quotient
    of multivariate polynomial ring in 2 variables (y//x), (z//x)
      over rational field
    by ideal (-(y//x)*(z//x) + 1)

julia> t = first(gens(OO(U)))
(t//s)

julia> Phi = morphism_from_rational_functions(IP1, C, U, V, [t//one(t), 1//t]);

julia> realizations = Oscar.realize_on_patch(Phi, U);

julia> realizations[3]
Affine scheme morphism
  from [(t//s)]          AA^1
  to   [(x//z), (y//z)]  scheme((x//z)^2 - (y//z))
given by the pullback function
  (x//z) -> (t//s)
  (y//z) -> (t//s)^2

```
"""
morphism_from_rational_functions(
      X::AbsCoveredScheme, Y::AbsCoveredScheme, 
      U::AbsAffineScheme, V::AbsAffineScheme,
      a::Vector{<:FieldElem};
      check::Bool=true,
      domain_covering::Covering=default_covering(X),
      codomain_covering::Covering=default_covering(Y)
     ) = MorphismFromRationalFunctions(X, Y, U, V, a; check, domain_covering, codomain_covering)

function Base.show(io::IOContext, Phi::MorphismFromRationalFunctions)
  if is_terse(io)
    print("Morphism from rational functions")
  else
    io = pretty(io)
    print(io, "Hom: ")
    print(io, Lowercase(), domain(Phi), " -> ", Lowercase(), codomain(Phi))
  end
end

function Base.show(io::IOContext, ::MIME"text/plain", Phi::MorphismFromRationalFunctions)
  io = pretty(io)
  println(io, "Morphism from rational functions")
  print(io, Indent())
  println(io, "from ", Lowercase(), domain(Phi))
  println(io, "to ", Lowercase(), codomain(Phi), Dedent())
  println(io, "with representatives")
  print(io, Indent())
  c = collect(patch_representatives(Phi))
  for (U,(V,imgs)) in c[1:end-1]
    print(io, "(")
    join(io, coordinates(V), ",")
    print(io, ") -> ")
    print(io, "(")
    join(io, imgs, ",")
    print(io, ")")
  end
  (U,(V,imgs)) = c[end]
  print(io, "(")
  join(io, coordinates(V), ",")
  print(io, ") -> ")
  print(io, "(")
  join(io, imgs, ",")
  print(io, ")")
  print(io, Dedent())
end
# For every pair of patches `U` in the `domain_covering` and `V` in the `codomain_covering` 
# the pullback `f₁,…,fᵣ` of the `gens` of `OO(V)` along `Phi` can be represented as 
# rational functions on `U`. This returns a dictionary where `U` can be used 
# as a key and a list of pairs `(V, [f₁,…,fᵣ])` is returned for every `V` for 
# which this has already been computed. 
patch_representatives(Phi::MorphismFromRationalFunctions) = Phi.patch_representatives

# The full realizations of the morphism: Keys are the `patches` `U` of the `domain_covering`
# and the output is a list of morphisms `φ : U' → V` from `PrincipalOpenSubset`s of `U` 
# to `patches` of the `codomain_covering` which are needed to provide a full 
# `CoveringMorphism` for `Phi`.
realizations(Phi::MorphismFromRationalFunctions) = Phi.realizations

# For every pair of patches `U` in the `domain_covering` and `V` in the `codomain_covering` 
# there is a maximal open subset `U' ⊂ U` (not necessarily principally open) so that 
# `φ : U' → V` is the restriction of `Phi` to `U'`. This returns a dictionary which takes 
# the pair `(U, V)` as input and returns a list of morphisms `φₖ : U'ₖ → V` with 
# all `U'ₖ` principally open in `U` and so that all the `U'ₖ` cover `U'`.
maximal_extensions(Phi::MorphismFromRationalFunctions) = Phi.maximal_extensions

# This is similar to `patch_representatives` only that this returns a dictionary 
# which takes pairs `(U, V)` as input and returns the pullback of `gens(OO(V))` as 
# rational functions in the fraction field of the `ambient_coordinate_ring` of `U`.
realization_previews(Phi::MorphismFromRationalFunctions) = Phi.realization_previews

# This is similar to `maximal_extensions`, but here only one `PrincipalOpenSubset` `U' ⊂ U` 
# is produced such that `Phi` can be realized as `φ : U' → V`, i.e. `U'` need not 
# be maximal with this property.
cheap_realizations(Phi::MorphismFromRationalFunctions) = Phi.cheap_realizations

@doc raw"""
    realize_on_patch(Phi::MorphismFromRationalFunctions, U::AbsAffineScheme)

For ``U`` in the `domain_covering` of `Phi` construct a list of morphisms 
``fₖ : U'ₖ → Vₖ`` from `PrincipalOpenSubset`s ``U'ₖ`` of ``U`` to `patches` 
``Vₖ`` in the `codomain_covering` so that altogether the `fₖ` can be assembled 
to a `CoveringMorphism` which realizes `Phi`.
"""
function realize_on_patch(Phi::MorphismFromRationalFunctions, U::AbsAffineScheme)
  if haskey(realizations(Phi), U)
    return realizations(Phi)[U]
  end
  X = domain(Phi)
  Y = codomain(Phi)
  V = codomain_chart(Phi)

  # Try to cover U by PrincipalOpenSubsets W so that the restriction 
  # of Phi to W extends to a regular morphism φ : W → V' for some 
  # `affine_chart` of the codomain of Phi.
  covered_codomain_patches = Vector{AbsAffineScheme}([V])
  complement_equations = Vector{elem_type(OO(U))}()
  FY = function_field(Y, check=Phi.run_internal_checks)
  FX = function_field(X, check=Phi.run_internal_checks)
  A = [FX(a) for a in coordinate_images(Phi)]
  a = [b[U] for b in A]
  #a = [lift(simplify(OO(U)(numerator(b))))//lift(simplify(OO(U)(denominator(b)))) for b in a]
  list_for_V = _extend(U, a)
  @assert !is_empty(list_for_V) "list must not be empty"
  Psi = [morphism(W, ambient_space(V), b, check=Phi.run_internal_checks) for (W, b) in list_for_V]
  # Up to now we have maps to the ambient space of V. 
  # But V might be a hypersurface complement in there and we 
  # might need to restrict our domain of definition accordingly. 
  Psi_res = AffineSchemeMor[_restrict_properly(psi, V; check=Phi.run_internal_checks) for psi in Psi]
  @assert all(phi->codomain(phi) === V, Psi_res)
  append!(complement_equations, [OO(U)(lifted_numerator(complement_equation(domain(psi)))) for psi in Psi_res])
  while !isone(ideal(OO(U), complement_equations))
    # Find another chart in the codomain which is hopefully easily accessible
    V_next, V_orig = _find_good_neighboring_patch(codomain_covering(Phi), covered_codomain_patches)
    # Get the gluing morphisms for the gluing to some already covered chart
    f, g = gluing_morphisms(gluings(codomain_covering(Phi))[(V_next, V_orig)])
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
    Psi = [morphism(W, ambient_space(V_next), b, check=Phi.run_internal_checks) for (W, b) in list_for_V_next]
    Psi = [_restrict_properly(psi, V_next; check=Phi.run_internal_checks) for psi in Psi]
    append!(Psi_res, Psi)
    append!(complement_equations, [OO(U)(lifted_numerator(complement_equation(domain(psi)))) for psi in Psi])
    push!(covered_codomain_patches, V_next)
  end
  realizations(Phi)[U] = Psi_res
  return Psi_res
end

@doc raw"""
    realize_on_open_subset(Phi::MorphismFromRationalFunctions, U::AbsAffineScheme, V::AbsAffineScheme)

Return a morphism `f : U' → V` from some `PrincipalOpenSubset` of `U` to `V` such
that the restriction of `Phi` to `U'` is `f`. Note that `U'` need not be maximal 
with this property!
"""
function realize_on_open_subset(Phi::MorphismFromRationalFunctions, U::AbsAffineScheme, V::AbsAffineScheme)
  X = domain(Phi)
  Y = codomain(Phi)
  # Check that the input is admissible
  if !any(x->x===U, patches(X)) 
    UU = _find_chart(U, default_covering(X))
  end
  if !any(x->x===V, patches(Y)) 
    VV = _find_chart(V, default_covering(Y))
  end
  y = function_field(Y; check=false).(gens(OO(V)))
  dom_rep = domain_chart(Phi)
  cod_rep = codomain_chart(Phi)
  y_cod = [a[cod_rep] for a in y]::Vector{<:FieldElem}
  x_dom = [evaluate(a, coordinate_images(Phi)) for a in y_cod]::Vector{<:FieldElem}
  x = function_field(X; check=false).(x_dom)
  img_gens_frac = [a[U] for a in x]
  dens = [denominator(a) for a in img_gens_frac]
  U_sub = PrincipalOpenSubset(U, OO(U).(dens))
  img_gens = [OO(U_sub)(numerator(a), denominator(a)) for a in img_gens_frac]
  prelim = morphism(U_sub, ambient_space(V), img_gens, check=Phi.run_internal_checks) # TODO: Set to false
  return _restrict_properly(prelim, V; check=Phi.run_internal_checks)
end

@doc raw"""
    realization_preview(Phi::MorphismFromRationalFunctions, U::AbsAffineScheme, V::AbsAffineScheme)

For a pair `(U, V)` of `patches` in the `domain_covering` and the `codomain_covering` 
of `Phi`, respectively, this returns a list of elements in the fraction field of the 
`ambient_coordinate_ring` of `U` which represent the pullbacks of `gens(OO(V))` under 
`Phi` to `U`.
"""
function realization_preview(Phi::MorphismFromRationalFunctions, U::AbsAffineScheme, V::AbsAffineScheme)
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
  y = function_field(Y; check=false).(gens(OO(V)))
  dom_rep = domain_chart(Phi)
  cod_rep = codomain_chart(Phi)
  y_cod = [a[cod_rep] for a in y]::Vector{<:FieldElem}
  x_dom = [evaluate(a, coordinate_images(Phi)) for a in y_cod]::Vector{<:FieldElem}
  x = function_field(X; check=false).(x_dom; check=true)
  img_gens_frac = [a[U] for a in x]
  realization_previews(Phi)[(U, V)] = img_gens_frac
  return img_gens_frac
end

@doc raw"""
    random_realization(Phi::MorphismFromRationalFunctions, U::AbsAffineScheme, V::AbsAffineScheme)

For a pair `(U, V)` of `patches` in the `domain_covering` and the `codomain_covering` 
of `Phi`, respectively, this creates a random `PrincipalOpenSubset` `U'` on which 
the restriction `f : U' → V` of `Phi` can be realized and returns that restriction.
Note that `U'` need not (and usually will not) be maximal with this property.
"""
function random_realization(Phi::MorphismFromRationalFunctions, U::AbsAffineScheme, V::AbsAffineScheme)
  img_gens_frac = realization_preview(Phi, U, V)
  U_sub, img_gens = _random_extension(U, img_gens_frac)
  phi = morphism(U_sub, ambient_space(V), img_gens, check=false)
  return phi
end

@doc raw"""
    cheap_realization(Phi::MorphismFromRationalFunctions, U::AbsAffineScheme, V::AbsAffineScheme)

For a pair `(U, V)` of `patches` in the `domain_covering` and the `codomain_covering` 
of `Phi`, respectively, this creates a random `PrincipalOpenSubset` `U'` on which 
the restriction `f : U' → V` of `Phi` can be realized and returns that restriction.
Note that `U'` need not (and usually will not) be maximal with this property.

This method is cheap in the sense that it simply inverts all representatives of 
the denominators occurring in the `realization_preview(Phi, U, V)`.
"""
function cheap_realization(Phi::MorphismFromRationalFunctions, U::AbsAffineScheme, V::AbsAffineScheme)
 #if haskey(cheap_realizations(Phi), (U, V))
 #  return cheap_realizations(Phi)[(U, V)]
 #end
  img_gens_frac = realization_preview(Phi, U, V)
  # Try to cancel the fractions heuristically; turns out it was too expensive in some applications due to slow divide
# for (k, f) in enumerate(img_gens_frac)
#   a = numerator(f)
#   b = denominator(f)
#   aa = OO(U)(a)
#   new_num = aa
#   new_den = one(aa)
#   fac = factor(b)
#   for (p, e) in fac
#     success, q = divides(aa, OO(U)(p))
#     while success && e > 0
#       aa = q
#       e = e - 1
#       success, q = divides(aa, OO(U)(p))
#     end
#     new_den = new_den * p^e
#   end
#   @assert aa*denominator(f) == unit(fac)*new_den*numerator(f)
#   img_gens_frac[k] = inv(unit(fac))*fraction(aa)//fraction(new_den)
# end
  denoms = OO(U).([denominator(a) for a in img_gens_frac])
  #any(is_zero, denoms) && error("some denominator was zero")
  U_sub = PrincipalOpenSubset(U, denoms)
  img_gens = [OO(U_sub)(numerator(a), denominator(a), check=false) for a in img_gens_frac] 
  phi = morphism(U_sub, ambient_space(V), img_gens, check=false) 
  cheap_realizations(Phi)[(U, V)] = phi
  return phi
end

@doc raw"""
    realize_maximally_on_open_subset(Phi::MorphismFromRationalFunctions, U::AbsAffineScheme, V::AbsAffineScheme)

For a pair `(U, V)` of `patches` in the `domain_covering` and the `codomain_covering` 
of `Phi`, respectively, this returns a list of morphisms `fₖ : U'ₖ → V` such that the 
restriction of `Phi` to `U'ₖ` and `V` is `fₖ` and altogether the `U'ₖ` cover the maximal 
open subset `U'⊂ U` on which the restriction `U' → V` of `Phi` can be realized.
"""
function realize_maximally_on_open_subset(Phi::MorphismFromRationalFunctions, U::AbsAffineScheme, V::AbsAffineScheme)
  if haskey(maximal_extensions(Phi), (U, V))
    return maximal_extensions(Phi)[(U, V)]
  end
  img_gens_frac = realization_preview(Phi, U, V)
  extensions = _extend(U, img_gens_frac)
  result = AbsAffineSchemeMor[]
  for (U, g) in extensions
    prelim = morphism(U, ambient_space(V), g, check=Phi.run_internal_checks)
    push!(result, _restrict_properly(prelim, V))
  end
  maximal_extensions(Phi)[(U, V)] = result
  return result
end


@doc raw"""
    realize(Phi::MorphismFromRationalFunctions)

Compute a full realization of `Phi` as a `CoveredSchemeMorphism`. Note 
that this computation is very expensive and usage of this method should 
be avoided.
"""
function realize(Phi::MorphismFromRationalFunctions; check::Bool=true)
  if !isdefined(Phi, :full_realization)
    realizations = AbsAffineSchemeMor[]
    mor_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
    for U in patches(domain_covering(Phi))
      loc_mors = realize_on_patch(Phi, U)
      for phi in loc_mors 
        mor_dict[domain(phi)] = phi
      end
      append!(realizations, loc_mors)
    end
    domain_ref = Covering([domain(phi) for phi in realizations])
    inherit_gluings!(domain_ref, domain_covering(Phi))
    # TODO: Inherit the decomposition_info, too!
    phi_cov = CoveringMorphism(domain_ref, codomain_covering(Phi), mor_dict; 
                               check=(Phi.run_internal_checks || check))
    # Make the refinement known to the domain
    push!(coverings(domain(Phi)), domain_ref)
    Phi.full_realization = CoveredSchemeMorphism(domain(Phi), codomain(Phi), phi_cov; 
                                                 check=(Phi.run_internal_checks || check))
  end
  return Phi.full_realization
end

underlying_morphism(Phi::MorphismFromRationalFunctions) = realize(Phi; check=Phi.run_internal_checks)

###
# Find a random open subset `W ⊂ U` to which all the rational functions 
# represented by the elements in `a` can be extended as regular functions.
function _random_extension(U::AbsAffineScheme, a::Vector{<:FieldElem})
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

###
# Find a maximal open subset `W ⊂ U` to which all the rational functions 
# represented by the elements in `a` can be extended as regular functions
# and return a list of tuples `(W', a')` of realizations on principal 
# open subsets W' covering W.
function _extend(U::AbsAffineScheme, a::Vector{<:FieldElem})
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

  result = Vector{Tuple{AbsAffineScheme, Vector{RingElem}}}()

  for g in small_generating_set(I_undef)
    Ug = PrincipalOpenSubset(U, g)
    b = [OO(Ug)(numerator(f), denominator(f)) for f in a]
    #b = [divides(OO(Ug)(numerator(f)), OO(Ug)(denominator(f)))[2] for f in a]
    push!(result, (Ug, b))
  end

  return result
end

# Some functionality that was missing and should probably be moved elsewhere.
# TODO: Do that.
equidimensional_decomposition_radical(I::MPolyLocalizedIdeal) = [ideal(base_ring(I), gens(J)) for J in equidimensional_decomposition_radical(saturated_ideal(I))]
equidimensional_decomposition_radical(I::MPolyQuoLocalizedIdeal) = [ideal(base_ring(I), gens(J)) for J in equidimensional_decomposition_radical(saturated_ideal(I))]

### When realizing a `MorphismFromRationalFunctions` `Phi` on pairs 
# of patches `(U, V)`, it is essential to use information on 
# other pairs `(U', V')` of patches which is already available 
# through feasible channels. Now, for example, for `U'` as above
# this finds another patch `U` for which the gluing of `U` and `U'` 
# is already fully computed, but which is not in `covered`. 
# If no such `U` exists: Bad luck. We just take any other one 
# and the gluing has to be computed eventually. 
function _find_good_neighboring_patch(cov::Covering, covered::Vector{<:AbsAffineScheme})
  U = [x for x in patches(cov) if !any(y->y===x, covered)]
  glue = gluings(cov)
  good_neighbors = [(x, y) for x in U for y in covered if 
                    haskey(glue, (x, y)) && 
                    (glue[(x, y)] isa SimpleGluing || 
                     (glue[(x, y)] isa LazyGluing && is_computed(glue[(x, y)]))
                    )
                   ]
  if !isempty(good_neighbors)
    return first(good_neighbors)
  end
  isempty(U) && error("no new neighbor could be found")
  return first(U), first(covered)
end

# Even though a list of rational functions might be realizable 
# as regular functions on U' and a morphism U' → A to the `ambient_space` 
# of V can be realized, V might be so small that we need a proper restriction 
# of the domain. The methods below take care of that. 
function _restrict_properly(f::AbsAffineSchemeMor, V::AbsAffineScheme{<:Ring, <:MPolyRing}; check::Bool=true)
  return restrict(f, domain(f), V; check)
end

function _restrict_properly(f::AbsAffineSchemeMor, V::AbsAffineScheme{<:Ring, <:MPolyQuoRing}; check::Bool=true)
  return restrict(f, domain(f), V; check)
end

function _restrict_properly(
    f::AbsAffineSchemeMor{<:PrincipalOpenSubset}, V::AbsAffineScheme{<:Ring, <:RT};
    check::Bool=true
  ) where {RT<:MPolyLocRing{<:Ring, <:RingElem, 
                            <:MPolyRing, <:MPolyRingElem, 
                            <:MPolyPowersOfElement}
          }
  h = denominators(inverted_set(OO(V)))
  pbh = pullback(f).(h)
  U = domain(f)
  W = ambient_scheme(U)
  UU = PrincipalOpenSubset(W, push!(OO(W).(lifted_numerator.(pbh)), complement_equation(U)))
  return restrict(f, UU, V; check)
end

function _restrict_properly(
    f::AbsAffineSchemeMor{<:PrincipalOpenSubset}, V::AbsAffineScheme{<:Ring, <:RT};
    check::Bool=true
  ) where {RT<:MPolyQuoLocRing{<:Ring, <:RingElem, 
                            <:MPolyRing, <:MPolyRingElem, 
                            <:MPolyPowersOfElement}
          }
  h = denominators(inverted_set(OO(V)))
  pbh = pullback(f).(h)
  U = domain(f)
  W = ambient_scheme(U)
  UU = PrincipalOpenSubset(W, push!(OO(W).(lifted_numerator.(pbh)), complement_equation(U)))
  return restrict(f, UU, V; check)
end

### The natural mathematical way to deal with algebraic cycles. However, since 
# we can not realize fraction fields of integral domains 𝕜[x₁,…,xₙ]/I properly, 
# not even to speak of their transcendence degrees, this functionality is rather 
# limited at the moment. 
function pushforward(Phi::MorphismFromRationalFunctions, D::AbsAlgebraicCycle)
  error("not implemented")
end

function pushforward(Phi::MorphismFromRationalFunctions, D::AbsWeilDivisor)
  is_isomorphism(Phi) || error("method not implemented unless for the case of an isomorphism")
  #is_proper(Phi) || error("morphism must be proper")
  all(is_prime, components(D)) || error("divisor must be given in terms of irreducible components")
  X = domain(Phi)
  Y = codomain(Phi)
  pushed_comps = IdDict{AbsIdealSheaf, elem_type(coefficient_ring(D))}()
  for I in components(D)
    J = _pushforward_prime_divisor(Phi, I) # Use dispatch here
    pushed_comps[J] = D[I]
  end
  is_empty(pushed_comps) && error("pushforward of this divisor along an alleged isomorphism is empty")
  return WeilDivisor(AlgebraicCycle(Y, coefficient_ring(D), pushed_comps); check=false)
end

# The following attributes can not be checked algorithmically at the moment. 
# But they can be set by the user so that certain checks of other methods 
# are satisfied; i.e. the user has to take responsibility and confirm that 
# they know what they're doing through these channels. 
@attr Bool function is_proper(phi::AbsCoveredSchemeMorphism)
  error("no method implemented to check properness")
end

@attr Bool function is_isomorphism(phi::AbsCoveredSchemeMorphism)
  error("no method implemented to check for being an isomorphism")
end
  
is_known(::typeof(is_isomorphism), phi::AbsCoveredSchemeMorphism) = has_attribute(phi, :is_isomorphism)

### Pullback of algebraic cycles along an isomorphism. 
function pullback(phi::MorphismFromRationalFunctions, C::AbsAlgebraicCycle)
  is_isomorphism(phi) || error("method is currently only implemented for isomorphisms")
  X = domain(phi)
  Y = codomain(phi)
  R = coefficient_ring(C)
  comps = IdDict{AbsIdealSheaf, elem_type(R)}()
  for I in components(C)
    @vprint :MorphismFromRationalFunctions 1 "trying cheap pullback\n"
    pbI = _try_pullback_cheap(phi, I)
    if pbI === nothing
      @vprint :MorphismFromRationalFunctions 1 "trying randomized pullback\n"
      pbI = _try_randomized_pullback(phi, I)
      if pbI === nothing
        @vprint :MorphismFromRationalFunctions 1 "trying the full pullback\n"
        pbI = _pullback(phi, I)
      end
    end
    comps[pbI] = C[I]
  end

  return AlgebraicCycle(X, R, comps)
end

# In order to pull back an ideal sheaf I along phi we need to find only pair of 
# dense open subsets (U, V) such that the restriction of `phi` can be realized 
# as a regular morphism f : U → V with f*(I) non-zero in OO(U). 
# The method below tries to find such a pair in a cheap way which might not 
# be successful.
function _try_pullback_cheap(phi::MorphismFromRationalFunctions, I::AbsIdealSheaf)
  X = domain(phi)
  Y = codomain(phi)
  scheme(I) === Y || error("ideal sheaf not defined on the correct scheme")
  # Find a patch in Y on which this component is visible
  all_V = [V for V in affine_charts(Y) if !isone(I(V))]
  function complexity_codomain(V::AbsAffineScheme)
    return sum(total_degree.(lifted_numerator.(gens(I(V)))); init=0)
  end
  sort!(all_V; by=complexity_codomain)
  for V in all_V

    # Find a patch in X in which the pullback is visible
    JJ = IdealSheaf(X)
    all_U = copy(affine_charts(X))
    function complexity(U::AbsAffineScheme)
      a = realization_preview(phi, U, V)
      return maximum(vcat([total_degree(numerator(f)) for f in a], [total_degree(denominator(f)) for f in a]))
    end
    sort!(all_U; by=complexity)

    # First try to get hold of the component via cheap realizations 
    pullbacks = IdDict{AbsAffineScheme, Ideal}()
    for U in all_U
      psi = cheap_realization(phi, U, V)
      U_sub = domain(psi)
      pullbacks[U] = pullback(psi)(saturated_ideal(I(V)))
    end
      #J = pullback(psi)(saturated_ideal(I(V)))
    function new_complexity(U::AbsAffineScheme)
      return sum(total_degree.(lifted_numerator.(gens(pullbacks[U]))); init=0)
    end
    sort!(all_U, lt=(x,y)->new_complexity(x)<new_complexity(y))
    for U in all_U
      J = pullbacks[U]
      psi = cheap_realization(phi, U, V)
      if !isone(J)
        JJ = IdealSheaf(X, domain(psi), gens(J))
        return JJ
        break
      end
    end
  end
  return nothing
end

#=
# The following was thought to be easier. But it turns out not to be.
# with a complicated map it is in general cheaper to first spread out 
# the ideal sheaf in the codomain and then have more choices as to which 
# pair of charts to use for pullback.
function _try_pullback_cheap(phi::MorphismFromRationalFunctions, I::PrimeIdealSheafFromChart)
  X = domain(phi)
  Y = codomain(phi)
  scheme(I) === Y || error("ideal sheaf not defined on the correct scheme")
  # Find a patch in Y on which this component is visible
  V0 = original_chart(I)

  for U in affine_charts(X)
    psi = cheap_realization(phi, U, V0)
    J = pullback(psi)(saturated_ideal(I(V0)))
    if !isone(J)
      JJ = IdealSheaf(X, domain(psi), gens(J))
      return JJ
      break
    end
  end

  return nothing
end

function _try_randomized_pullback(phi::MorphismFromRationalFunctions, I::PrimeIdealSheafFromChart)
  X = domain(phi)
  Y = codomain(phi)
  scheme(I) === Y || error("ideal sheaf not defined on the correct scheme")
  # Find a patch in Y on which this component is visible
  V0 = original_chart(I)

  for U in affine_charts(X)
    psi = random_realization(phi, U, V0)
    J = pullback(psi)(saturated_ideal(I(V0)))
    if !isone(J)
      JJ = IdealSheaf(X, domain(psi), gens(J))
      return JJ
      break
    end
  end

  return nothing
end
=#

function _pullback(phi::MorphismFromRationalFunctions, I::PrimeIdealSheafFromChart)
  V0 = original_chart(I)
  X = domain(phi)
  Y = codomain(phi)
  V = affine_charts(Y)
  U = affine_charts(X)

  cod_patches = AbsAffineScheme[V0]
  cod_patches = vcat(cod_patches, [U for U in keys(object_cache(I)) if any(x->x===U, affine_charts(Y))])
  cod_patches = vcat(cod_patches, [U for U in affine_charts(Y) if !any(x->x===U, cod_patches)])
  for U0 in U
    I_undef = ideal(OO(U0), elem_type(OO(U0))[])
    random_realizations = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
    for V0 in cod_patches
      psi = random_realization(phi, U0, V0)
      I_undef = I_undef + ideal(OO(U0), complement_equation(domain(psi)))
      random_realizations[V0] = psi
      if isone(I_undef)
        break
      end
    end
    if isone(I_undef)
      J = ideal(OO(U0), elem_type(OO(U0))[])
      for (V0, psi) in random_realizations
        J_loc = pullback(psi, saturated_ideal(I(V0)))
        J = J + ideal(OO(U0), lifted_numerator.(gens(J_loc)))
      end
      return IdealSheaf(X, U0, gens(J))
    else
      continue
    end
  end
  error("pullback did not succeed")
end

# Similar to the above function, but this time we try pairs (U, V) and determine the 
# maximal open subset W ⊂ U such that the restriction `W → V` of `phi` can be realized. 
# Then we take a random linear combination `h` of the generators of the ideal for the 
# complement of W in U and realize the restriction of `phi` on the hypersurface complement 
# of `h`. With probability 1 this will produce a non-trivial pullback of I on this 
# patch whenever I was non-trivial on V. But it is not as cheap as the method above 
# since the rational functions must be converted to regular functions on D(h). 
function _try_randomized_pullback(phi::MorphismFromRationalFunctions, I::AbsIdealSheaf)
  X = domain(phi)
  Y = codomain(phi)
  scheme(I) === Y || error("ideal sheaf not defined on the correct scheme")
  # Find a patch in Y on which this component is visible
  all_V = [V for V in affine_charts(Y) if !isone(I(V))]

  min_var = minimum([ngens(OO(V)) for V in all_V])
  all_V = [V for V in all_V if ngens(OO(V)) == min_var]
  deg_bound = minimum([maximum([total_degree(lifted_numerator(g)) for g in gens(I(V))]) for V in all_V])
  all_V = [V for V in all_V if maximum([total_degree(lifted_numerator(g)) for g in gens(I(V))]) == deg_bound]
  V = first(all_V)

  all_U = copy(affine_charts(X))
  function complexity(U::AbsAffineScheme)
    a = realization_preview(phi, U, V)
    return maximum(vcat([total_degree(numerator(f)) for f in a], [total_degree(denominator(f)) for f in a]))
  end
  sort!(all_U, by=complexity)

  for U in all_U
    psi = random_realization(phi, U, V)

    J = pullback(psi)(saturated_ideal(I(V)))
    if !isone(J)
      JJ = IdealSheaf(X, domain(psi), gens(J))
      return JJ
    end
  end
  return nothing
end

### Deprecated method below, left here for recycling.
function _pullback(phi::MorphismFromRationalFunctions, I::AbsIdealSheaf)
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
  function complexity(U::AbsAffineScheme)
    a = realization_preview(phi, U, V)
    return maximum(vcat([total_degree(numerator(f)) for f in a], [total_degree(denominator(f)) for f in a]))
  end
  sort!(all_U; by=complexity)

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

function pullback(phi::MorphismFromRationalFunctions, D::AbsWeilDivisor)
  return WeilDivisor(pullback(phi)(underlying_cycle(D)), check=false) 
end

function _find_good_representative_chart(I::PrimeIdealSheafFromChart)
  return original_chart(I)
end

function _find_good_representative_chart(I::PullbackIdealSheaf)
  f = morphism(I)
  f_cov = covering_morphism(f)
  J = original_ideal_sheaf(I)
  V = _find_good_representative_chart(J)
  list = maps_with_given_codomain(f_cov, V)
  for f_loc in list
    !isone(I(domain(f_loc))) && return domain(f_loc)
  end

  # if the above doesn't help, fall back to the default
  X = scheme(I)
  for U in keys(object_cache(I))
    any(x->x===U, affine_charts(X)) || continue
    !is_one(I(U)) && return U
  end
  for U in affine_charts(X)
    !is_one(I(U)) && return U
  end
  error("no chart found")
end

function _find_good_representative_chart(I::AbsIdealSheaf; covering::Covering=default_covering(scheme(I)))
  # We assume that I is prime
  # TODO: Make this an hassert?
  @hassert :IdealSheaves 2 is_prime(I)
  X = scheme(I)

  # Some heuristics to choose a reasonably "easy" chart
  cand = AbsAffineScheme[]
  for U in keys(object_cache(I))
    any(x->x===U, patches(covering)) || continue
    !is_one(I(U)) && push!(cand, U)
  end

  function complexity(U::AbsAffineScheme)
    g = lifted_numerator.(gens(I(U)))
    return maximum(total_degree(f) for f in g; init=0)
  end

  if !is_empty(cand)
    c = complexity.(cand)
    m = minimum(c)
    i = findfirst(==(m), c)
    return cand[i]
  end

  for U in patches(covering)
    !is_one(I(U)) && return U
  end
  error("no chart found")
end

function _prepare_pushforward_prime_divisor(
    phi::MorphismFromRationalFunctions, I::AbsIdealSheaf;
    domain_chart::AbsAffineScheme = _find_good_representative_chart(I),
    codomain_charts::Vector{<:AbsAffineScheme} = copy(patches(codomain_covering(phi)))
  )
  @assert !is_one(I(domain_chart))
  U = domain_chart
  X = domain(phi)
  Y = codomain(phi)

  # try cheap realizations first
  sorted_charts = copy(codomain_charts)
  if has_decomposition_info(default_covering(Y))
    info = decomposition_info(default_covering(Y))
    # Enabling the following line seems to lead to wrong results. Why?
    #sorted_charts = filter!(V->dim(OO(V)) - dim(ideal(OO(V), elem_type(OO(V))[OO(V)(a) for a in info[V]])) <= 1, sorted_charts)
  end

  function compl(V::AbsAffineScheme)
    result = 0
    if (U, V) in keys(realization_previews(phi))
      fracs = realization_previews(phi)[(U, V)]::Vector
      if any(f->OO(U)(denominator(f)) in I(U), fracs)
        result = result + 100000
      else
        #result = sum(length(terms(numerator(f))) + length(terms(denominator(f))) for f in fracs; init=0)
        result = sum(total_degree(numerator(f)) + total_degree(denominator(f)) for f in fracs; init=0)
      end
    else
      result = result + 10
    end
    return result
  end

  sort!(sorted_charts; by=compl)

  bad_charts = Int[]
  for (i, V) in enumerate(sorted_charts)
    # Find a chart in the codomain which has a chance to have the pushforward visible
    fracs = realization_preview(phi, U, V)::Vector
    any(f->OO(U)(denominator(f)) in I(U), fracs) && continue
    phi_loc = cheap_realization(phi, U, V)
    # Shortcut to decide whether the restriction will lead to a trivial ideal
    if OO(V) isa MPolyLocRing || OO(V) isa MPolyQuoLocRing
      is_bad_chart = false
      for h in denominators(inverted_set(OO(V)))
        if pullback(phi_loc)(h) in I(domain(phi_loc)) 
          # Remove this chart from the list
          push!(bad_charts, i)
          is_bad_chart = true
          break
        end
      end
      is_bad_chart && continue
    end
    return phi_loc, U, V
  end

  sorted_charts = AbsAffineScheme[V for (i, V) in enumerate(sorted_charts) if !(i in bad_charts)]
  sort!(sorted_charts; by=compl)
  
  # try random realizations second
  loc_ring, _ = localization(OO(U), complement_of_prime_ideal(I(U)))

  # The ring is smooth in codimension one. Let's find a generator of its maximal ideal
  pp = ideal(loc_ring, gens(I(U)))
  qq = pp^2
  candidates = [g for g in gens(I(U)) if !(loc_ring(g) in qq)]
  complexity(a) = total_degree(lifted_numerator(a)) + total_degree(lifted_denominator(a))
  sort!(candidates, by=complexity)
  isempty(candidates) && error("no element of valuation one found")

  min_terms = minimum(length.(terms.(lifted_numerator.(candidates))))
  h = candidates[findfirst(x->length(terms(lifted_numerator(x)))==min_terms, candidates)]

  F1 = FreeMod(loc_ring, 1)

  # Trigger caching of the attribute :is_prime for faster computation of is_zero
  # on elements.
  if loc_ring isa MPolyQuoLocRing
    is_prime(modulus(underlying_quotient(loc_ring)))
  end

  P, _ = sub(F1, [h*F1[1]]) # The maximal ideal in the localized ring, but as a submodule

  for V in sorted_charts
    fs = realization_preview(phi, U, V)
    skip = false
    for (i, fr) in enumerate(fs)
      a = numerator(fr)
      b = denominator(fr)
      aa = loc_ring(a)
      bb = loc_ring(b)
      count = 0
      # If the denominator is in P, we have a problem.
      # If the numerator is not in P, the problem is serious and this chart can 
      # not be used.
      # If the numerator is also in P, we can cancel the fraction by h and 
      # start all over. 
      while OO(U)(lifted_numerator(bb)) in I(U)
        count = count + 1
        if !(OO(U)(lifted_numerator(aa)) in I(U))
          skip = true
          break
        end
        bb = coordinates(bb*F1[1], P)[1]
        aa = coordinates(aa*F1[1], P)[1]
      end
      skip && break
      num = lifted_numerator(aa)*lifted_denominator(bb)
      den = lifted_numerator(bb)*lifted_denominator(aa)
      @assert !(OO(U)(den) in I(U))
      fs[i] = num//den
    end
    skip && continue

    # Copied from cheap_realization
    denoms = OO(U).([denominator(a) for a in fs])
    any(is_zero, denoms) && error("some denominator was zero in this chart")
    U_sub = PrincipalOpenSubset(U, denoms)
    img_gens = [OO(U_sub)(numerator(a), denominator(a), check=false) for a in fs]
    psi = morphism(U_sub, ambient_space(V), img_gens, check=false)
    # TODO: Do we really want to cache this? The expressions become more complex by the above cancellation.
    #cheap_realizations(phi)[(U, V)] = psi
    #realization_previews(phi)[(U, V)] = fs

    # Shortcut to decide whether the restriction will lead to a trivial ideal
    if OO(V) isa MPolyLocRing || OO(V) isa MPolyQuoLocRing
      is_bad_chart = false
      I_sub = ideal(OO(U_sub), lifted_numerator.(gens(I(U)))) # Avoid production of the ring map, etc.
      for h in denominators(inverted_set(OO(V)))
        if OO(U_sub)(pullback(psi)(h)) in I_sub
          is_bad_chart = true
          break
        end
      end
      is_bad_chart && continue
    end

    return psi, U, V
    # Else: try the next chart
  end

  # The preselection of charts in the codomain via the optional argument 
  # may lead to that there is no result in the end. 
  return nothing, Oscar.domain_chart(phi), codomain_chart(phi)
end

function _pushforward_prime_divisor(
    phi::MorphismFromRationalFunctions, D::AbsWeilDivisor
  )
  R = coefficient_ring(D)
  div_dict = IdDict{AbsIdealSheaf, elem_type(R)}()
  for I in components(D)
    div_dict[_pushforward_prime_divisor(phi, I)] = D[I]
  end
  return WeilDivisor(AlgebraicCycle(domain(phi), R, div_dict; check=false); check=false)
end


function _pushforward_prime_divisor(
    phi::MorphismFromRationalFunctions, I::AbsIdealSheaf;
    domain_chart::AbsAffineScheme=_find_good_representative_chart(I),
    codomain_charts::Vector{<:AbsAffineScheme} = copy(patches(codomain_covering(phi)))
  )
  loc_map, dom_chart, cod_chart = _prepare_pushforward_prime_divisor(phi, I; domain_chart, codomain_charts)
  loc_map === nothing && return nothing
  U_sub = domain(loc_map)
  J = preimage(pullback(loc_map), I(U_sub))
  JJ = ideal(OO(cod_chart), gens(J))
  return PrimeIdealSheafFromChart(codomain(phi), cod_chart, JJ)
end


function compose(
    f::MorphismFromRationalFunctions,
    g::MorphismFromRationalFunctions
  )
  @assert codomain(f) === domain(g)
  fracs = coordinate_images(g)
  X = domain(f)
  Y = domain(g)
  Z = codomain(g)
  FY = function_field(Y; check=false)
  imgs = FY.(fracs)
  V = codomain_chart(f)
  imgs_V = [b[V] for b in imgs]
  U = domain_chart(f)
  imgs_U = [evaluate(numerator(h), coordinate_images(f))//evaluate(denominator(h), coordinate_images(f)) for h in imgs_V]
  fg =  morphism_from_rational_functions(X, Z, U, codomain_chart(g), imgs_U; check=false)
  if is_known(is_isomorphism, f) && is_known(is_isomorphism, g)
    if is_isomorphism(f) && is_isomorphism(g)
      set_attribute!(fg, :is_isomorphism=>true)
    end 
  end
  return fg 
end

#=
# This will not work because `phi` must be dominant. But then its image 
# has the zero ideal sheaf associated to it.
function ideal_sheaf_of_image(phi::MorphismFromRationalFunctions)
  X = domain(phi)
  Y = codomain(phi)
  U = domain_chart(phi)
  V = codomain_chart(phi)
  charts = vcat([V], affine_charts(Y))
  for V in charts
    phi_loc = random_realization(phi, U, V)
    pb_phi = pullback(phi_loc)
    K = kernel(pb_phi)
    res = ideal(OO(V), elem_type(OO(V))[OO(V)(a) for a in gens(K)])
    if !is_one(res) # TODO: Abbreviate test via complement equations?
      return PrimeIdealSheafFromChart(Y, V, res)
    end
  end
  error("image ideal could not be computed")
end
=#

