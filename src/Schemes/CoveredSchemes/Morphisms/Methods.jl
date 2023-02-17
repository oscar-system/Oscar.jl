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
  set_attribute!(X, :simplified_covering, Csimp)
  return X
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
  if codomain(cf) === domain(cg) 
    mor_dict = IdDict{AbsSpec, AbsSpecMor}() # TODO: Keep the type of the morphisms?
    for U in basic_patches(domain(cf))
      mor_dict[U] = compose(cf[U], cg[codomain(cf[U])])
    end
    cc = CoveringMorphism(domain(cf), codomain(cg), mor_dict, check=false)
    return CoveredSchemeMorphism(X, Z, cc, check=false)
  else
    haskey(refinements(Y), (codomain(cf), domain(cg))) || error("composition of this complicated case is not yet implemented")
    ref_mor = Y[codomain(cf), domain(cg)] # the refinement `CoveringMorphism`
    mor_dict = IdDict{AbsSpec, AbsSpecMor}() # TODO: Keep the type of the morphisms?
    for U in basic_patches(domain(cf))
      V = codomain(cf[U])
      W = codomain(ref_mor[V])
      mor_dict[U] = compose(compose(cf[U], ref_mor[V]), cg[W])
    end
    cc = CoveringMorphism(domain(cf), codomain(cg), mor_dict, check=false)
    return CoveredSchemeMorphism(X, Z, cc, check=false)
  end
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

