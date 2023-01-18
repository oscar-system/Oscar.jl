@Markdown.doc """
    pullback(f::CoveredSchemeMorphism)

Return a function `phi` which takes any reasonable 
argument `v` associated to the `codomain` of `f` and 
produces ``f*(v)``.
"""
function pullback(f::CoveredSchemeMorphism)
  X = domain(f)
  Y = codomain(f)

  function pullback_func(II::IdealSheaf)
    scheme(II) === Y || error("ideal sheaf is not defined on the codomain of the function")
    phi = covering_morphism(f)
    ID = IdDict{AbsSpec, Ideal}()
    for U in patches(domain(phi))
      f_U = phi[U]
      V = codomain(f_U)
      pbf = pullback(f_U)
      ID[U] = ideal(OO(U), pbf.(gens(II(V))))
    end
    return IdealSheaf(X, ID, check=false) # TODO: Set to false eventually.
  end

  function pullback_func(C::CartierDivisor)
    phi = covering_morphism(f)
    triv_dict = IdDict{AbsSpec, RingElem}()
    OOX = OO(X)
    # We need to do the following:
    # The divisor D := f^* C needs to be trivialized on some refinement. 
    # The one we're using consists of the preimages of the patches of 
    # the `trivializing_covering` of `C`. These must be intersected 
    # with each and every affine chart of X first. 
    for U in patches(domain(phi))
      V = codomain(phi[U])
      for W in patches(trivializing_covering(C))
        # In case V and W do not happen to lay in the same affine 
        # chart, the preimage of W in U will be empty. So we can 
        # skip the computation in that case.
        success, par = _have_common_ancestor(V, W)
        if success
          # Reconstruct W as a PrincipalOpenSubset of the affine chart directly
          # (Note that we might have gone through some calls to `simplify(...)`, 
          # so this really is a non-trivial step)
          inc_W_flat = _flatten_open_subscheme(W, par)
          # To cheaply construct the preimage of W in U, just pull back the 
          # complement equation.
          h = complement_equation(codomain(inc_W_flat))
          UW = PrincipalOpenSubset(U, pullback(phi[U])(OOX(par, V)(h)))
          # Finally, we can add the trivialization on U ∩ f⁻¹(W) to our list.
          triv_dict[UW] = pullback(phi[U])(first(C(V)))
        end
      end
    end
    return CartierDivisor(X, triv_dict) # We have to leave the checks back on, because local 
                                        # equations might pull back to zero or zero divisors.
  end

  return pullback_func
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

  return SimpleGlueing(X, Y, ff, gg)
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

