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

