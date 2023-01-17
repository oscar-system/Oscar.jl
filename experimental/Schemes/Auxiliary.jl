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
    for U in patches(domain(phi))
      V = codomain(phi[U])
      for W in patches(trivializing_covering(C))
        success, par = _have_common_ancestor(V, W)
        if success
          inc_W_flat = _flatten_open_subscheme(W, par)
          h = complement_equation(codomain(inc_W_flat))
          UW = PrincipalOpenSubset(U, pullback(phi[U])(OOX(par, V)(h)))
          triv_dict[UW] = pullback(phi[U])(first(C(V)))
        end
      end
    end
    return CartierDivisor(X, triv_dict) # We have to leave the checks back on, because local 
                                        # equations might pull back to zero or zero divisors.
  end

  return pullback_func
end

