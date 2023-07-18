########################################################################
# Constructors for CoveringMorphism                                    #
########################################################################

@doc raw"""
    identity_map(C::Covering) -> CoveringMorphism

Given a covering `C`, return the covering morphism given as the identity
on each patch of `C`.

# Examples
```jldoctest
julia> P, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> Y = variety(ideal([x^3-y^2*z]));

julia> Ycov = covered_scheme(Y);

julia> C = default_covering(Ycov)
Covering
  described by patches
    1: spec of quotient of multivariate polynomial ring
    2: spec of quotient of multivariate polynomial ring
    3: spec of quotient of multivariate polynomial ring
  in the coordinate(s)
    1: [(y//x), (z//x)]
    2: [(x//y), (z//y)]
    3: [(x//z), (y//z)]

julia> identity_map(C)
Morphism
  from covering with 3 patches
    1a: [(y//x), (z//x)]   spec of quotient of multivariate polynomial ring
    2a: [(x//y), (z//y)]   spec of quotient of multivariate polynomial ring
    3a: [(x//z), (y//z)]   spec of quotient of multivariate polynomial ring
  to   covering with 3 patches
    1b: [(y//x), (z//x)]   spec of quotient of multivariate polynomial ring
    2b: [(x//y), (z//y)]   spec of quotient of multivariate polynomial ring
    3b: [(x//z), (y//z)]   spec of quotient of multivariate polynomial ring
given by the pullback functions
  1a -> 1b
    (y//x) -> (y//x)
    (z//x) -> (z//x)
  ----------------------------------------
  2a -> 2b
    (x//y) -> (x//y)
    (z//y) -> (z//y)
  ----------------------------------------
  3a -> 3b
    (x//z) -> (x//z)
    (y//z) -> (y//z)
```
"""
function identity_map(C::Covering)
  map_dict = IdDict{AbsSpec, AbsSpecMor}()
  for U in patches(C)
    map_dict[U] = identity_map(U)
  end
  return CoveringMorphism(C, C, map_dict, check=false)
end


