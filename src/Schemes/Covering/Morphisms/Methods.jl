########################################################################
# Composition of CoveringMorphisms                                     #
########################################################################
function compose(f::CoveringMorphism, g::CoveringMorphism)
  domain(g) === codomain(f) || error("morphisms can not be composed")
  morphism_dict = IdDict{AbsSpec, AbsSpecMor}()
  for U in patches(domain(f))
    morphism_dict[U] = compose(f[U], g[codomain(f[U])])
  end
  return CoveringMorphism(domain(f), codomain(g), morphism_dict)
end

########################################################################
# Simplification of Coverings                                          #
########################################################################
@Markdown.doc """
    simplify(C::Covering)

Given a covering ``C`` apply `simplify` to all basic affine patches 
in ``C`` and return a triple ``(C', f, g)`` consisting of the 
resulting covering ``C'`` and the identifying isomorphism 
``f : C' ↔ C``.
"""
function simplify(C::Covering)
  n = npatches(C)
  new_patches = [simplify(X) for X in patches(C)]
  GD = glueings(C)
  new_glueings = IdDict{Tuple{AbsSpec, AbsSpec}, AbsGlueing}()
  for (X, Y) in keys(GD)
    Xsimp = new_patches[C[X]]
    iX, jX = identification_maps(Xsimp)
    Ysimp = new_patches[C[Y]]
    iY, jY = identification_maps(Ysimp)
    G = GD[(X, Y)]
    new_glueings[(Xsimp, Ysimp)] = restrict(G, jX, jY, check=false)
  end
  iDict = IdDict{AbsSpec, AbsSpecMor}()
  jDict = IdDict{AbsSpec, AbsSpecMor}()
  for i in 1:length(new_patches)
    iDict[new_patches[i]] = identification_maps(new_patches[i])[1]
    jDict[C[i]] = identification_maps(new_patches[i])[2]
  end
  Cnew = Covering([ U for U in new_patches], new_glueings, check=false)
  i_cov_mor = CoveringMorphism(Cnew, C, iDict)
  j_cov_mor = CoveringMorphism(C, Cnew, jDict)
  return Cnew, i_cov_mor, j_cov_mor
end

