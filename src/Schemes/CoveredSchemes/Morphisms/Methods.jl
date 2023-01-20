export simplify

########################################################################
# Methods for CoveredSchemeMorphism                                    #
########################################################################
@Markdown.doc """
    simplify!(X::AbsCoveredScheme)

Apply `simplify` to the `default_covering` of `X` and store the 
resulting `Covering` ``C'`` and its identification with the default 
covering in `X`; return a tuple ``(X, C')``.
"""
function simplify!(X::AbsCoveredScheme)
  C = default_covering(X)
  Csimp, i, j = simplify(C)
  push!(coverings(X), Csimp)
  refinements(X)[(C, Csimp)] = j
  refinements(X)[(Csimp, C)] = i
  return (X, Csimp)
end

########################################################################
# Auxiliary methods for compatibility                                  #
########################################################################
lifted_numerator(f::MPolyElem) = f
lifted_numerator(f::MPolyQuoElem) = lift(f)

function compose(f::AbsCoveredSchemeMorphism, g::AbsCoveredSchemeMorphism)
  X = domain(f)
  Y = codomain(f)
  Y === domain(g) || error("maps not compatible")
  Z = codomain(g)
  cf = covering_morphism(f)
  cg = covering_morphism(g)
  codomain(cf) === domain(cg) || error("change of coverings not implemented yet")
  mor_dict = IdDict{AbsSpec, AbsSpecMor}() # TODO: Keep the type of the morphisms?
  for U in basic_patches(domain(cf))
    mor_dict[U] = compose(cf[U], cg[codomain(cf[U])])
  end
  cc = CoveringMorphism(domain(cf), codomain(cg), mor_dict, check=false)
  return CoveredSchemeMorphism(X, Z, cc, check=false)
end

function maps_with_given_codomain(f::AbsCoveredSchemeMorphism, V::AbsSpec)
  fcov = covering_morphism(f)
  result = Vector{AbsSpecMor}()
  for U in keys(morphisms(fcov))
    floc = morphisms(fcov)[U]
    codomain(floc) === V || continue
    push!(result, floc)
  end
  return result
end

