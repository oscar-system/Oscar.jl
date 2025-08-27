########################################################################
# production and restriction maps for the structure sheaf
########################################################################

### objects
function produce_object(F::StructureSheafOfRings, U::AbsAffineScheme)
  return OO(U)
end

function produce_object(F::StructureSheafOfRings, U::AffineSchemeOpenSubscheme)
  return OO(U)
end

### restrictions
function produce_restriction_map(F::StructureSheafOfRings, V::AbsAffineScheme, U::AbsAffineScheme)
  V === U || error("basic affine patches must be the same")
  return id_hom(OO(V))
end

function produce_restriction_map(
    F::StructureSheafOfRings, 
    V::AbsAffineScheme, 
    U::Union{<:PrincipalOpenSubset, <:SimplifiedAffineScheme}
  )
  X = scheme(F)
  OV = F(V) # Assumed to be cached or produced on the fly.
  OU = F(U) # Same as above.
  incU = _flatten_open_subscheme(U, default_covering(X))
  #incU, dU = _find_chart(U, default_covering(X))
  U_flat = codomain(incU)
  W = ambient_scheme(U_flat)
  if W === V
    return pullback(compose(incU, inclusion_morphism(U_flat)))
  else
    G = default_covering(X)[V, W]
    f, g = gluing_morphisms(G)
    pbg = pullback(g)
    function rho_func(x::RingElem)
      parent(x) === OV || error("element does not belong to the correct domain")
      y = pbg(domain(pbg)(x, check=false))
      yy = restrict(y, U_flat, check=false)
      return pullback(incU)(yy)
    end
    return hom(OV, OU, rho_func.(gens(OV)), check=false)
  end
  error("arguments are not valid")
end

function produce_restriction_map(F::StructureSheafOfRings, 
    V::Union{<:PrincipalOpenSubset, <:SimplifiedAffineScheme},
    U::AbsAffineScheme
  )
  X = scheme(F)
  OV = F(V)
  OU = F(U) 
  incV = _flatten_open_subscheme(V, default_covering(X))
  W = ambient_scheme(codomain(incV))
  V_direct = domain(incV)
  if W === U
    # By virtue of the checks in _is_open_func we must have V isomorphic to U.
    phi = pullback(inverse(incV))
    psi = hom(OO(V_direct), OU, gens(OU))
    return hom(OV, OU, psi.(phi.(gens(OV))))
    ### deprecated code below;
    # kept for the moment because of possible incompatibilities with gluings 
    # along AffineSchemeOpenSubschemes.
    function rho_func(a::RingElem)
      parent(a) === OV || error("element does not belong to the correct ring")
      # We may assume that all denominators admissible in V are
      # already units in OO(U)
      return OU(lifted_numerator(a))*inv(OU(lifted_denominator(a)))
    end
    return hom(OV, OU, rho_func.(gens(OV)), check=false)
  else
    G = default_covering(X)[W, U]
    W1, W2 = gluing_domains(G)
    f, g = gluing_morphisms(G)
    g_res = restrict(g, U, V_direct, check=false)
    return pullback(compose(g_res, inverse(incV)))
    ### deprecated code below; see comment above
    function rho_func2(a::RingElem)
      parent(a) === OV || error("element does not belong to the correct ring")
      return restrict(pullback(g)(OO(W1)(a)), U)
    end
    return hom(OV, OU, rho_func2.(gens(OV)), check=false)
  end
end
function produce_restriction_map(F::StructureSheafOfRings, 
    V::Union{<:PrincipalOpenSubset, <:SimplifiedAffineScheme},
    U::Union{<:PrincipalOpenSubset, <:SimplifiedAffineScheme}
  )
  X = scheme(F)
  OV = F(V)
  OU = F(U)
  inc_U_flat = _flatten_open_subscheme(U, default_covering(X))
  inc_V_flat = _flatten_open_subscheme(V, default_covering(X))
  A = ambient_scheme(codomain(inc_U_flat))
  B = ambient_scheme(codomain(inc_V_flat))
  U_flat = codomain(inc_U_flat)
  V_flat = codomain(inc_V_flat)

  if A === B
    return hom(OV, OU, 
               pullback(inc_U_flat).(pullback(inclusion_morphism(U_flat, V_flat, check=false)).(pullback(inverse(inc_V_flat)).(gens(OV)))), check=false
              )
  else
    G = default_covering(X)[A, B]
    f, g = gluing_morphisms(G)
    VV_flat = intersect(V_flat, codomain(f))
    VU = preimage(f, VV_flat, check=false)
    fres = restrict(f, VU, VV_flat, check=false)
    inc_V_flat_inv = inverse(inc_V_flat)
    function rho_func(x::RingElem)
      parent(x) === OV || error("input not valid")
      y = pullback(inverse(inc_V_flat))(x)
      y = restrict(y, VV_flat, check=false)
      y = pullback(fres)(y)
      y = restrict(y, U_flat, check=false)
      return pullback(inc_U_flat)(y)
    end
    return hom(OV, OU, rho_func.(gens(OV)), check=false)
  end
  error("arguments are invalid")
end

function produce_restriction_map(F::StructureSheafOfRings, V::AbsAffineScheme, W::AffineSchemeOpenSubscheme)
  X = scheme(F)
  OV = F(V)
  OW = F(W)
  V in default_covering(X) || return false
  ambient_scheme(W) in default_covering(X) || return false
  if V === ambient_scheme(W)
    return MapFromFunc(OV, OW, OW)
  else
    G = default_covering(X)[V, ambient_scheme(W)]
    f, g = gluing_morphisms(G)
    function rho_func(a::RingElem)
      parent(a) === OV || error("element does not belong to the correct ring")
      return restrict(pullback(g)(OO(domain(f))(a)), W, check=false)
    end
    return MapFromFunc(OV, OW, rho_func)
  end
end

### cleaned up until here ###
# We do not make AffineSchemeOpenSubscheme compatible with the tree structures, yet. 
# All AffineSchemeOpenSubscheme's are hence required to have an ambient_scheme on the top level. 

function produce_restriction_map(F::StructureSheafOfRings, V::PrincipalOpenSubset, W::AffineSchemeOpenSubscheme)
  error("method not implemented at the moment")
  X = scheme(F)
  OV = F(V)
  OW = F(W)
  if ambient_scheme(V) === ambient_scheme(W)
    function rho_func(a::RingElem)
      parent(a) === OV || error("element does not belong to the correct ring")
      return OW(a)
    end
    return MapFromFunc(OV, OW, rho_func)
  else
    G = default_covering(X)(ambient_scheme(V), ambient_scheme(W))
    f, g = gluing_morphisms(G)
    VG = intersect(V, domain(f))
    preV = preimage(g, VG, check=false)
    gres = restriction(g, preV, VG, check=false)
    inc = inclusion_morphism(W, preV, check=false)
    function rho_func2(a::RingElem)
      parent(a) === OV || error("element does not belong to the correct ring")
      return pullback(inc)(pullback(gres)(OO(preV)(a)))
    end
    return MapFromFunc(OV, OW, rho_func2)
  end
end

function produce_restriction_map(F::StructureSheafOfRings, V::AffineSchemeOpenSubscheme, W::AffineSchemeOpenSubscheme)
  X = scheme(F)
  OV = F(V)
  OW = F(W)
  if ambient_scheme(V) === ambient_scheme(W)
    inc = inclusion_morphism(W, V, check=false)
    return MapFromFunc(OV, OW, pullback(inc))
  else
    G = default_covering(X)[ambient_scheme(V), ambient_scheme(W)]
    f, g = gluing_morphisms(G)
    VG = intersect(V, domain(f))
    inc0 = inclusion_morphism(VG, V, check=false)
    preV = preimage(g, VG, check=false)
    gres = restrict(g, preV, VG, check=false)
    inc = inclusion_morphism(W, preV, check=false)
    return MapFromFunc(OV, OW, x->(pullback(inc)(pullback(gres)(pullback(inc0)(x)))))
  end
end

