export covering_morphism, domain_covering, codomain_covering

########################################################################
# Interface for AbsCoveredSchemeMorphism                               #
########################################################################
### essential getters 
function domain(f::AbsCoveredSchemeMorphism)::CoveredScheme
  return domain(underlying_morphism(f))
end

function codomain(f::AbsCoveredSchemeMorphism)::CoveredScheme
  return codomain(underlying_morphism(f))
end

function covering_morphism(f::AbsCoveredSchemeMorphism)::CoveringMorphism
  return covering_morphism(underlying_morphism(f))
end

### generically derived getters
domain_covering(f::AbsCoveredSchemeMorphism) = domain(covering_morphism(f))
codomain_covering(f::AbsCoveredSchemeMorphism) = codomain(covering_morphism(f))
getindex(f::AbsCoveredSchemeMorphism, U::Spec) = covering_morphism(f)[U]

########################################################################
# Basic getters for CoveredSchemeMorphism                              #
########################################################################
domain(f::CoveredSchemeMorphism) = f.X
codomain(f::CoveredSchemeMorphism) = f.Y
covering_morphism(f::CoveredSchemeMorphism) = f.f

