
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

@doc raw"""
    isomorphism_on_open_subsets(f::AbsCoveredSchemeMorphism)

For a birational morphism ``f : X → Y`` of `AbsCoveredScheme`s this 
returns an isomorphism of affine schemes ``fᵣₑₛ : U → V`` which is 
the restriction of ``f`` to two dense open subsets ``U ⊂ X`` and 
``V ⊂ Y``.
"""
function isomorphism_on_open_subsets(f::AbsCoveredSchemeMorphism)
  if !has_attribute(f, :iso_on_open_subset)
    is_birational(f) # Should compute and store the attribute
  end
  return get_attribute(f, :iso_on_open_subset)::AbsSpecMor
end

