########################################################################
# Composition of CoveringMorphisms                                     #
########################################################################
function compose(f::CoveringMorphism, g::CoveringMorphism)
  domain(g) === codomain(f) || error("morphisms can not be composed")
  morphism_dict = IdDict{AbsSpec, AbsSpecMor}()
  for U in patches(domain(f))
    morphism_dict[U] = compose(f[U], g[codomain(f[U])])
  end
  return CoveringMorphism(domain(f), codomain(g), morphism_dict, check=false)
end

########################################################################
# Simplification of Coverings                                          #
########################################################################
@doc raw"""
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
    #new_glueings[(Xsimp, Ysimp)] = restrict(G, jX, jY, check=false)
    new_glueings[(Xsimp, Ysimp)] = LazyGlueing(Xsimp, Ysimp, _compute_restriction, 
                                               RestrictionDataIsomorphism(G, jX, jY)
                                              )
  end
  iDict = IdDict{AbsSpec, AbsSpecMor}()
  jDict = IdDict{AbsSpec, AbsSpecMor}()
  for i in 1:length(new_patches)
    iDict[new_patches[i]] = identification_maps(new_patches[i])[1]
    jDict[C[i]] = identification_maps(new_patches[i])[2]
  end
  Cnew = Covering([ U for U in new_patches], new_glueings, check=false)
  i_cov_mor = CoveringMorphism(Cnew, C, iDict, check=false)
  j_cov_mor = CoveringMorphism(C, Cnew, jDict, check=false)
  return Cnew, i_cov_mor, j_cov_mor
end

########################################################################
# Base change
########################################################################
function base_change(phi::Any, f::CoveringMorphism;
    domain_map::CoveringMorphism=base_change(phi, domain(f))[2],
    codomain_map::CoveringMorphism=base_change(phi, codomain(f))[2]
  )
  D = domain(f)
  C = codomain(f)
  DD = domain(domain_map)
  CC = domain(codomain_map)
  mor_dict = IdDict{AbsSpec, AbsSpecMor}()
  for UU in patches(DD)
    U = codomain(domain_map[UU])
    V = codomain(f[U])
    g_V = first(maps_with_given_codomain(codomain_map, V)) # The result must be unique as it arises 
                                                           # from a base change.
    _, ff, _ = base_change(phi, f[U], domain_map=domain_map[UU], codomain_map=g_V)
    mor_dict[UU] = ff
  end

  return domain_map, CoveringMorphism(DD, CC, mor_dict, check=true), codomain_map # TODO: Set to false after testing.
end

