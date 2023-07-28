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
defined by ideal(y, x)

julia> Xcov = covered_scheme(X)
Scheme
  over rational field
with default covering
  described by patches
    1: spec of quotient of multivariate polynomial ring
    2: spec of quotient of multivariate polynomial ring
    3: spec of quotient of multivariate polynomial ring
  in the coordinate(s)
    1: [(y//x), (z//x)]
    2: [(x//y), (z//y)]
    3: [(x//z), (y//z)]

julia> identity_map(Xcov)
Morphism
  from scheme over QQ covered with 3 patches
    1a: [(y//x), (z//x)]   spec of quotient of multivariate polynomial ring
    2a: [(x//y), (z//y)]   spec of quotient of multivariate polynomial ring
    3a: [(x//z), (y//z)]   spec of quotient of multivariate polynomial ring
  to   scheme over QQ covered with 3 patches
    1b: [(y//x), (z//x)]   spec of quotient of multivariate polynomial ring
    2b: [(x//y), (z//y)]   spec of quotient of multivariate polynomial ring
    3b: [(x//z), (y//z)]   spec of quotient of multivariate polynomial ring
given by the pullback functions
  1a -> 1b
    (y//x) -> 0
    (z//x) -> 0
    ----------------------------------------
  2a -> 2b
    (x//y) -> 0
    (z//y) -> 0
    ----------------------------------------
  3a -> 3b
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
