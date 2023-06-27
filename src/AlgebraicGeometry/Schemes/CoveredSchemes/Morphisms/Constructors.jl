########################################################################
# Constructors for CoveredSchemeMorphism                               #
########################################################################

### This type has no external constructors.
@attr AbsCoveredSchemeMorphism function identity_map(X::AbsCoveredScheme)
  C = default_covering(X)
  id_cov = identity_map(C)
  result = CoveredSchemeMorphism(X, X, id_cov)
  set_attribute!(result, :inverse, result)
  return result
end
