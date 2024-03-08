########################################################################
# Constructors for CoveredSchemeMorphism                               #
########################################################################

### This type has no external constructors.

@doc raw"""
    identity_map(X::AbsCoveredScheme) -> AbsCoveredSchemeMorphism

Given a covered scheme `X`, return the identity map of `X` seen as a
covered scheme morphism.

# Examples
```jldoctest
julia> P, x = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> X = variety(ideal(x[1:2]))
Projective variety
  in projective 2-space over QQ with coordinates [x, y, z]
defined by ideal (y, x)

julia> Xcov = covered_scheme(X)
Scheme
  over rational field
with default covering
  described by patches
    1: scheme((y//z), (x//z))
  in the coordinate(s)
    1: [(x//z), (y//z)]

julia> identity_map(Xcov)
Covered scheme morphism
  from scheme over QQ covered with 1 patch
    1a: [(x//z), (y//z)]   scheme((y//z), (x//z))
  to scheme over QQ covered with 1 patch
    1b: [(x//z), (y//z)]   scheme((y//z), (x//z))
given by the pullback function
  1a -> 1b
    (x//z) -> 0
    (y//z) -> 0
```
"""
@attr AbsCoveredSchemeMorphism function identity_map(X::AbsCoveredScheme)
  C = default_covering(X)
  id_cov = identity_map(C)
  result = CoveredSchemeMorphism(X, X, id_cov)
  set_attribute!(result, :inverse, result)
  return result
end
