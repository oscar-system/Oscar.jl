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
end

