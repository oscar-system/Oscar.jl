# This struct associates to a morphism f of type `MorphismType` 
# the mathematical symbol f^* which, in one way or the other, 
# denotes some pullback functor. 
#
# Which functor exactly this should be, is determined only from 
# the context in which this symbol appears. Hence, it is no 
# rigorous mathematical object, but really only a symbol. 
#
# It is introduced to facilitate writing 
#
#   pbf = pullback(f) # f some morphism f : X → Y
#   N = pbf(M) # M some object on Y
# 
# rather than 
# 
#   N = pullback(f, M)
#
# In particular, this makes the pullback a univariate 
# function which might come in handy when passing it as an argument 
# to other functions. For instance, for a collection of objects 
# A = [M₁, M₂, … ] on Y you can now write 
#
#   B = pbf.(A)
#
# rather than
#
#   B = (a->pullback(f, a).(A)
#
# or
#
#   B = [pullback(f, a) for a in A]
struct UniversalPullbackSymbol{MorphismType}
  f::MorphismType
end

# The following refers a call of f_star(M) to a bivariate pullback-method.
# The latter is what really needs to be implemented by the programmer. 
# Also, all parent checks, etc. are expected to happen there.
function (f_star::UniversalPullbackSymbol)(M::Any)
  return pullback(f_star.f, M)
end

@Markdown.doc """
    pullback(f::CoveredSchemeMorphism)

Return a function `phi` which takes any reasonable 
argument `M` associated to the `codomain` of `f` and 
produces ``f*(M)``.

**Note:** Internally, this simply calls `pullback(f, M)`. 
Hence, that method needs to be implemented.
"""
function pullback(f::AbsCoveredSchemeMorphism)
  return UniversalPullbackSymbol{typeof(f)}(f)
end

function pullback(f::AbsCoveredSchemeMorphism, II::IdealSheaf)
  X = domain(f)
  Y = codomain(f)
  scheme(II) === Y || error("ideal sheaf is not defined on the codomain of the function")
  phi = covering_morphism(f)
  ID = IdDict{AbsSpec, Ideal}()
  for U in patches(domain(phi))
    f_U = phi[U]
    V = codomain(f_U)
    pbf = pullback(f_U)
    ID[U] = ideal(OO(U), pbf.(gens(II(V))))
  end
  return IdealSheaf(X, ID, check=false)
end

function pullback(f::AbsCoveredSchemeMorphism, C::EffectiveCartierDivisor)
  X = domain(f)
  Y = codomain(f)
  phi = covering_morphism(f)
  triv_dict = IdDict{AbsSpec, RingElem}()
  psi = restrict(f, trivializing_covering(C))
  triv_cov = domain(psi)
  for U in patches(triv_cov)
    V = codomain(psi[U])
    pbgens = pullback(psi[U]).(C(V))

    # Do the sanity checks
    length(pbgens) == 1 || error("cartier divisor is not principal on this patch")
    g = first(pbgens)
    # See the Stacks project on cartier divisors
    is_zero_divisor(g) && error("pullback of the local equation is a zero divisor; pullback of the cartier divisor is not defined")
    triv_dict[U] = g
  end

  return EffectiveCartierDivisor(X, triv_dict, trivializing_covering=triv_cov, check=false)
end

function pullback(f::AbsCoveredSchemeMorphism, CC::Covering)
  psi = restrict(f, CC)
  return domain(psi)
end

function pullback(f::AbsCoveredSchemeMorphism, M::AbsCoherentSheaf)
  return PullbackSheaf(f, M)
end

@Markdown.doc """
    restrict(f::CoveredSchemeMorphism, DD::Covering)

Given a morphism of `CoveredScheme`s ``f : X → Y``, described via ``φ : C → D`` 
on `Covering`s `C` of ``X`` and `D` of ``Y``, and a refinement `DD ≤ D`, 
compute the preimages ``Uᵢ ∩ f⁻¹(Vⱼ)`` for ``Uᵢ`` the `patches` in `C` and 
``Vⱼ`` those of `DD` and restrictions of ``f`` to these new patches. 
Return the resulting `CoveringMorphism` ``ψ : C ∩ f⁻¹(DD) → DD``.
"""
function restrict(f::AbsCoveredSchemeMorphism, DD::Covering)
  X = domain(f)
  Y = codomain(f)
  phi = covering_morphism(f)
  C = domain(phi)
  D = codomain(phi)
  # DD needs to be a refinement of D; otherwise quit.
  all(x->some_ancestor(y->any(z->(z===y), patches(D)), x), patches(DD)) || error("second argument needs to be a refinement of the codomain covering on which the first argument is defined")

  res_cache = restriction_cache(f)
  haskey(res_cache, DD) && return res_cache[DD]

  OOX = OO(X)
  # We need to do the following:
  # The covering CC has patches Vⱼ in the codomain Y of f.
  # Their preimages must be taken in every patch Uᵢ of X in 
  # the domain's covering for phi. 
  res_dict = IdDict{AbsSpec, AbsSpecMor}()
  for U in patches(domain(phi))
    V = codomain(phi[U])
    for W in patches(DD)
      # In case V and W do not happen to lay in the same affine 
      # chart, the preimage of W in U will be empty. So we can 
      # skip the computation in that case.
      success, par = _have_common_ancestor(V, W)
      if success
        # Reconstruct W as a PrincipalOpenSubset of the affine chart directly
        # (Note that we might have gone through some calls to `simplify(...)`, 
        # so this really is a non-trivial step)
        iso_W_flat = _flatten_open_subscheme(W, par)
        W_flat = codomain(iso_W_flat)
        # To cheaply construct the preimage of W in U, just pull back the 
        # complement equation.
        h = complement_equation(codomain(iso_W_flat))
        UW = PrincipalOpenSubset(U, pullback(phi[U])(OOX(par, V)(h)))
        # Manually assemble the restriction of phi to this new patch
        ff = SpecMor(UW, W_flat, 
                     hom(OO(W_flat), OO(UW), OO(UW).(pullback(phi[U]).(gens(OO(V)))), check=false),
                     check=false
                    )
        f_res = compose(ff, inverse(iso_W_flat))
        res_dict[UW] = f_res
      end
    end
  end
  new_domain = Covering(collect(keys(res_dict)), IdDict{Tuple{AbsSpec, AbsSpec}, AbsGlueing}())
  inherit_glueings!(new_domain, domain(phi))

  psi = CoveringMorphism(new_domain, DD, res_dict, check=true)
  restriction_cache(f)[DD] = psi
  return psi
end

# Several objects on the codomain of f : X → Y might use the same covering 
# in their internals. Now when pulling them back, we do not want to recreate 
# the necessary refinement on X again and again. In particular, because 
# we do not want the pulled back modules/ideals/functions... to live on different 
# affine patches in X while their originals were defined on the same patches in Y.
@attr IdDict{<:Covering, <:CoveringMorphism} function restriction_cache(f::CoveredSchemeMorphism)
  res_dict = IdDict{Covering, CoveringMorphism}()
  phi = covering_morphism(f)
  res_dict[codomain(phi)] = phi
  return res_dict
end

struct InheritGlueingData
  orig::Covering
  X::AbsSpec
  Y::AbsSpec
end

function _compute_inherited_glueing(gd::InheritGlueingData)
  X = gd.X
  Y = gd.Y
  C = gd.orig

  success, Z = _have_common_ancestor(X, Y)
  if success
    # This is the easy case: Glueing within one chart. 
    # Keep in mind, however, that we might have gone through some simplify(...) calls, 
    # so we can not assume everything to be happening in the same ambient_ring.
    iso_X = _flatten_open_subscheme(X, Z)
    iso_Y = _flatten_open_subscheme(Y, Z)
    h_X = complement_equation(codomain(iso_X))
    h_Y = complement_equation(codomain(iso_Y))
    XY = PrincipalOpenSubset(X, pullback(iso_X)(OO(codomain(iso_X))(h_Y)))
    YX = PrincipalOpenSubset(Y, pullback(iso_Y)(OO(codomain(iso_Y))(h_X)))
    XYZ = PrincipalOpenSubset(Z, h_X*h_Y)

    x_img = gens(OO(X))
    x_img = pullback(inverse(iso_X)).(x_img)
    x_img = OO(XYZ).(x_img)
    phi = restrict(iso_Y, YX, XYZ)
    x_img = pullback(phi).(x_img)
    g = SpecMor(YX, XY, hom(OO(XY), OO(YX), x_img))
    
    y_img = gens(OO(Y))
    y_img = pullback(inverse(iso_Y)).(y_img)
    y_img = OO(XYZ).(y_img)
    psi = restrict(iso_X, XY, XYZ)
    y_img = pullback(psi).(y_img)
    f = SpecMor(XY, YX, hom(OO(YX), OO(XY), y_img))
    return SimpleGlueing(X, Y, f, g)
  end

  # As the easy case would have been caught before, we are now facing an inherited 
  # glueing across charts: X ↪ A ⊃ U ≅ V ⊂ B ↩ Y. 
  # We need to compute the intersection of X with Y along the identifications of U and V
  # and cook up the SimpleGlueing from that.
  iso_X = _flatten_open_subscheme(X, C)
  iso_Y = _flatten_open_subscheme(Y, C)
  A = ambient_scheme(codomain(iso_X))
  B = ambient_scheme(codomain(iso_Y))
  G = C[A, B] # The original glueing needed
  U, V = glueing_domains(G)
  f, g = glueing_morphisms(G)
  U isa PrincipalOpenSubset && ambient_scheme(U) === A || error("incorrect intermediate result")
  V isa PrincipalOpenSubset && ambient_scheme(V) === B || error("incorrect intermediate result")

  UX = intersect(U, codomain(iso_X))
  UX isa PrincipalOpenSubset && ambient_scheme(UX) === A || error("incorrect intermediate result")
  h_Y = pullback(f)(complement_equation(codomain(iso_Y)))
  UXY = PrincipalOpenSubset(codomain(iso_X), 
                            OO(codomain(iso_X))(lifted_numerator(h_Y)*complement_equation(U)))
  UXY isa PrincipalOpenSubset && ambient_scheme(UXY) === codomain(iso_X) || error("incorrect intermediate output")
  UY = PrincipalOpenSubset(U, h_Y)
  XY = PrincipalOpenSubset(X, pullback(iso_X)(complement_equation(UXY)))

  VY = intersect(V, codomain(iso_Y))
  VY isa PrincipalOpenSubset && ambient_scheme(VY) === B || error("incorrect intermediate result")
  h_X = pullback(g)(complement_equation(codomain(iso_X)))
  VX = PrincipalOpenSubset(V, h_X)
  VYX = PrincipalOpenSubset(codomain(iso_Y), 
                            OO(codomain(iso_Y))(lifted_numerator(h_X)*complement_equation(V)))
  VYX isa PrincipalOpenSubset && ambient_scheme(VYX) === codomain(iso_Y) || error("incorrect intermediate output")
  YX = PrincipalOpenSubset(Y, pullback(iso_Y)(complement_equation(VYX)))

  fres = restrict(f, UXY, VYX)
  gres = restrict(g, VYX, UXY)

  x_img = gens(OO(X))
  x_img = pullback(inverse(iso_X)).(x_img)
  x_img = OO(UXY).(x_img)
  x_img = pullback(gres).(x_img)
  phi = restrict(iso_Y, YX, VYX)
  x_img = pullback(phi).(x_img)
  gg = SpecMor(YX, XY, hom(OO(XY), OO(YX), x_img))

  y_img = gens(OO(Y))
  y_img = pullback(inverse(iso_Y)).(y_img)
  y_img = OO(VYX).(y_img)
  y_img = pullback(fres).(y_img)
  psi = restrict(iso_X, XY, UXY)
  y_img = pullback(psi).(y_img)
  ff = SpecMor(XY, YX, hom(OO(YX), OO(XY), y_img))

  return SimpleGlueing(X, Y, ff, gg, check=false)
end

@Markdown.doc """
    inherit_glueings!(ref::Covering, orig::Covering)

For a refinement `ref` of a `Covering` `orig` add all missing glueings 
in `ref` as inherited glueings from those in `orig`.
"""
function inherit_glueings!(ref::Covering, orig::Covering)
  for U in patches(ref)
    for V in patches(ref)
      if !haskey(glueings(ref), (U, V))
        glueings(ref)[(U, V)] = LazyGlueing(U, V, 
                                            _compute_inherited_glueing,
                                            InheritGlueingData(orig, U, V)
                                           )
      end
    end
  end
  return ref
end

