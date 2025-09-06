########################################################################
# Constructors for CoveringMorphism                                    #
########################################################################

@doc raw"""
    id_hom(C::Covering) -> CoveringMorphism

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
    1: scheme(-(y//x)^2*(z//x) + 1)
    2: scheme((x//y)^3 - (z//y))
    3: scheme((x//z)^3 - (y//z)^2)
  in the coordinate(s)
    1: [(y//x), (z//x)]
    2: [(x//y), (z//y)]
    3: [(x//z), (y//z)]

julia> id_hom(C)
Covering morphism
  from covering with 3 patches
    1a: [(y//x), (z//x)]   scheme(-(y//x)^2*(z//x) + 1)
    2a: [(x//y), (z//y)]   scheme((x//y)^3 - (z//y))
    3a: [(x//z), (y//z)]   scheme((x//z)^3 - (y//z)^2)
  to covering with 3 patches
    1b: [(y//x), (z//x)]   scheme(-(y//x)^2*(z//x) + 1)
    2b: [(x//y), (z//y)]   scheme((x//y)^3 - (z//y))
    3b: [(x//z), (y//z)]   scheme((x//z)^3 - (y//z)^2)
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
function id_hom(C::Covering)
  map_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
  for U in patches(C)
    map_dict[U] = id_hom(U)
  end
  return CoveringMorphism(C, C, map_dict, check=false)
end

@doc raw"""
    refinement_morphism(ref::Covering, orig::Covering)

Suppose `ref` is a refinement of `orig` just by implicit ancestry 
of the patches (every patch `U` of `ref` is a `PrincipalOpenSubset` 
of some `PrincipalOpenSubset` of some ... of some patch `V` in orig).

This function produces the inclusion map `ref -> orig` which is 
realized on all `patches` as `PrincipalOpenEmbedding`s. 
"""
function refinement_morphism(ref::Covering, orig::Covering)
  map_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
  for U in patches(ref)
    inc, h = _find_chart(U, orig)
    # TODO: construct and cache the inverse on image
    map_dict[U] = PrincipalOpenEmbedding(inc, h; check=false)
  end
  return CoveringMorphism(ref, orig, map_dict; check=false)
end
