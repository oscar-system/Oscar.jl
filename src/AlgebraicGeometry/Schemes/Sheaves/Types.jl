########################################################################
# Sheaves                                                              #
########################################################################
@doc raw"""
    AbsPreSheaf{SpaceType, OpenType, OutputType, RestrictionType}

Abstract type for a sheaf ℱ on a space X.

 * `SpaceType` is a parameter for the type of the space ``X`` on which ``ℱ`` is defined.

 * `OpenType` is a type (most probably abstract!) for the open sets ``U ⊂ X`` which are admissible as input for ``ℱ(U)``.

 * `OutputType` is a type (most probably abstract!) for the values that ``ℱ`` takes on admissible open sets ``U``.

 * `RestrictionType` is a parameter for the type of the restriction maps ``ℱ(V) → ℱ(U)`` for ``U ⊂ V ⊂ X`` open.

For any instance `F` of `AbsPreSheaf` on a topological space `X` the following methods are implemented:

 * `F(U)` for *admissible* open subsets ``U ⊂ X``: This returns the value ``ℱ(U)`` of the sheaf `F` on `U`. Note that due to technical limitations, not every type of open subset might be admissible.

 * `restriction_map(F, U, V)` for *admissible* open subsets ``V ⊂ U ⊂ X``: This returns the restriction map ``ρ : ℱ(U) → ℱ(V)``. Alternatively, one may also call `F(U, V)` to get this map.
"""
abstract type AbsPreSheaf{SpaceType, OpenType, OutputType, RestrictionType} end

########################################################################
# A minimal implementation of the sheaf interface on a scheme          #
########################################################################

@doc raw"""
    PreSheafOnScheme

A basic minimal implementation of the interface for `AbsPreSheaf`; to be used internally.
"""
@attributes mutable struct PreSheafOnScheme{SpaceType, OpenType, OutputType, RestrictionType,
                                      } <: AbsPreSheaf{
                                       SpaceType, OpenType,
                                       OutputType, RestrictionType
                                      }
  X::SpaceType

  # caches
  obj_cache::IdDict{<:OpenType, <:OutputType} # To cache values that have already been computed
  #res_cache::IdDict{<:Tuple{<:OpenType, <:OpenType}, <:RestrictionType} # To cache already computed restrictions

  is_open_func::Function # To check whether one set is open in the other

  function PreSheafOnScheme(X::Scheme; 
      OpenType=AbsAffineScheme, OutputType=Any, RestrictionType=Any,
      is_open_func::Any=is_open_embedding
    )
    return new{typeof(X), OpenType, OutputType, RestrictionType,
              }(X, IdDict{OpenType, OutputType}(),
                #IdDict{Tuple{OpenType, OpenType}, RestrictionType}(),
                is_open_func
               )
  end
end

########################################################################
# The structure sheaf of affine and covered schemes                    #
########################################################################
@doc raw"""
    StructureSheafOfRings <: AbsPreSheaf

On an `AbsCoveredScheme` ``X`` this returns the sheaf ``𝒪`` of rings of
regular functions on ``X``.

Note that due to technical reasons, the admissible open subsets are restricted
to the following:
 * `U::AbsAffineScheme` among the `basic_patches` of the `default_covering` of `X`;
 * `U::PrincipalOpenSubset` with `ambient_scheme(U)` in the `basic_patches` of the `default_covering` of `X`;
 * `W::AffineSchemeOpenSubscheme` with `ambient_scheme(W)` in the `basic_patches` of the `default_covering` of `X`.

One can call the restriction maps of ``𝒪`` across charts, implicitly using the
identifications given by the gluings in the `default_covering`.
# Examples
```jldoctest
julia> IP2 = projective_space(GF(7), [:x, :y, :z])
Projective space of dimension 2
  over prime field of characteristic 7
with homogeneous coordinates [x, y, z]

julia> X = covered_scheme(IP2)
Scheme
  over prime field of characteristic 7
with default covering
  described by patches
    1: affine 2-space
    2: affine 2-space
    3: affine 2-space
  in the coordinate(s)
    1: [(y//x), (z//x)]
    2: [(x//y), (z//y)]
    3: [(x//z), (y//z)]

julia> OOX = OO(X)
Structure sheaf of rings of regular functions
  on scheme over GF(7) covered with 3 patches
    1: [(y//x), (z//x)]   affine 2-space
    2: [(x//y), (z//y)]   affine 2-space
    3: [(x//z), (y//z)]   affine 2-space

julia> typeof(OOX)
StructureSheafOfRings{CoveredScheme{FqField}, Union{AbsAffineScheme, AffineSchemeOpenSubscheme}, Ring, Map}

julia> U, V, W = affine_charts(X)
3-element Vector{AffineScheme{FqField, FqMPolyRing}}:
 Affine 2-space
 Affine 2-space
 Affine 2-space

julia> glue = default_covering(X)[U, V]
Gluing
  of affine 2-space
  and affine 2-space
along the open subsets
  [(y//x), (z//x)]   AA^2 \ scheme((y//x))
  [(x//y), (z//y)]   AA^2 \ scheme((x//y))
given by the pullback function
  (x//y) -> 1/(y//x)
  (z//y) -> (z//x)/(y//x)

julia> UV, VU = gluing_domains(glue)
(AA^2 \ scheme((y//x)), AA^2 \ scheme((x//y)))

julia> y, z = gens(OOX(U))
2-element Vector{FqMPolyRingElem}:
 (y//x)
 (z//x)

julia> pb = OOX(U, VU)
Ring homomorphism
  from multivariate polynomial ring in 2 variables over GF(7)
  to localization of multivariate polynomial ring in 2 variables over GF(7) at products of ((x//y))
defined by
  (y//x) -> 1/(x//y)
  (z//x) -> (z//y)/(x//y)

julia> pb(y^2)
1/(x//y)^2

```
"""
@attributes mutable struct StructureSheafOfRings{SpaceType, OpenType, OutputType,
                                                 RestrictionType
                                                } <: AbsPreSheaf{
                                                                 SpaceType, OpenType,
                                                                 OutputType, RestrictionType
                                                                }
  OO::PreSheafOnScheme

  ### Structure sheaf on affine schemes
  function StructureSheafOfRings(X::AbsAffineScheme)
    function is_open_func(U::AbsAffineScheme, V::AbsAffineScheme)
      return is_subset(V, X) && is_open_embedding(U, V) # Note the restriction to subsets of X
    end
    R = PreSheafOnScheme(X,
                    OpenType=AbsAffineScheme, OutputType=Ring,
                    RestrictionType=Map,
                    is_open_func=is_open_func
                   )
    return new{typeof(X), Union{AbsAffineScheme, AffineSchemeOpenSubscheme}, Ring, Map}(R)
  end

  ### Structure sheaf on covered schemes
  function StructureSheafOfRings(X::AbsCoveredScheme)
    R = PreSheafOnScheme(X,
                      OpenType=Union{AbsAffineScheme, AffineSchemeOpenSubscheme}, OutputType=Ring,
                      RestrictionType=Map,
                      is_open_func=_is_open_func_for_schemes(X)
                     )
    return new{typeof(X), Union{AbsAffineScheme, AffineSchemeOpenSubscheme}, Ring, Map}(R)
  end
end

########################################################################
# Ideal sheaves on covered schemes                                     #
########################################################################

@doc raw"""
    AbsIdealSheaf <: AbsPreSheaf

A sheaf of ideals ``I`` on an `AbsCoveredScheme` ``X``.

For an affine open subset ``U ⊂ X`` call ``I(U)`` to obtain an ideal 
in `OO(U)` representing `I`.
# Examples
```jldoctest
julia> IP2 = projective_space(GF(7), [:x, :y, :z])
Projective space of dimension 2
  over prime field of characteristic 7
with homogeneous coordinates [x, y, z]

julia> X = covered_scheme(IP2)
Scheme
  over prime field of characteristic 7
with default covering
  described by patches
    1: affine 2-space
    2: affine 2-space
    3: affine 2-space
  in the coordinate(s)
    1: [(y//x), (z//x)]
    2: [(x//y), (z//y)]
    3: [(x//z), (y//z)]

julia> S = homogeneous_coordinate_ring(IP2)
Multivariate polynomial ring in 3 variables over GF(7) graded by
  x -> [1]
  y -> [1]
  z -> [1]

julia> II = ideal_sheaf(IP2, ideal(S, [S[1] + S[2]]))
Sheaf of ideals
  on scheme over GF(7) covered with 3 patches
    1: [(y//x), (z//x)]   affine 2-space
    2: [(x//y), (z//y)]   affine 2-space
    3: [(x//z), (y//z)]   affine 2-space
with restrictions
  1: Ideal ((y//x) + 1)
  2: Ideal ((x//y) + 1)
  3: Ideal ((x//z) + (y//z))

julia> U, V, W = affine_charts(X)
3-element Vector{AffineScheme{FqField, FqMPolyRing}}:
 Affine 2-space
 Affine 2-space
 Affine 2-space

julia> II(U)
Ideal generated by
  (y//x) + 1

julia> II(V)
Ideal generated by
  (x//y) + 1

julia> glue = default_covering(X)[U, V];

julia> UV, VU = gluing_domains(glue);

julia> II(U, VU) # transition functions are ring homomorphisms for ideal sheaves!
Ring homomorphism
  from multivariate polynomial ring in 2 variables over GF(7)
  to localization of multivariate polynomial ring in 2 variables over GF(7) at products of ((x//y))
defined by
  (y//x) -> 1/(x//y)
  (z//x) -> (z//y)/(x//y)

```
"""
abstract type AbsIdealSheaf{SpaceType, OpenType, OutputType,
                            RestrictionType
                           } <: AbsPreSheaf{
                                            SpaceType, OpenType,
                                            OutputType, RestrictionType
                                           }
end

@doc raw"""
    IdealSheaf <: AbsIdealSheaf

A sheaf of ideals ``ℐ`` on an `AbsCoveredScheme` ``X`` which is specified 
by a collection of concrete ideals on some open covering of ``X``.
"""
@attributes mutable struct IdealSheaf{SpaceType, OpenType, OutputType,
                                      RestrictionType
                                     } <: AbsIdealSheaf{
                                                      SpaceType, OpenType,
                                                      OutputType, RestrictionType
                                                     }
  ID::IdDict{AbsAffineScheme, Ideal} # the ideals on the patches of some covering of X
  OOX::StructureSheafOfRings # the structure sheaf on X
  I::PreSheafOnScheme # the underlying presheaf of ideals for caching

  ### Ideal sheaves on covered schemes
  function IdealSheaf(X::AbsCoveredScheme, ID::IdDict{AbsAffineScheme, Ideal};
      check::Bool=true
    )
    Ipre = PreSheafOnScheme(X, 
                      OpenType=AbsAffineScheme, OutputType=Ideal,
                      RestrictionType=Map,
                      is_open_func=_is_open_func_for_schemes_without_affine_scheme_open_subscheme(X)
                     )
    I = new{typeof(X), AbsAffineScheme, Ideal, Map}(ID, OO(X), Ipre)
    @check begin
      # Check that all ideal sheaves are compatible on the overlaps.
      # TODO: eventually replace by a check that on every basic
      # affine patch, the ideal sheaf can be inferred from what is
      # given on one dense open subset.
      C = default_covering(X)
      for U in basic_patches(default_covering(X))
        for V in basic_patches(default_covering(X))
          G = C[U, V]
          A, B = gluing_domains(G)
          for i in 1:number_of_complement_equations(A)
            I(A[i]) == ideal(OO(X)(A[i]), I(V, A[i]).(gens(I(V)))) || error("ideals do not coincide on overlap")
          end
          for i in 1:number_of_complement_equations(B)
            I(B[i]) == ideal(OO(X)(B[i]), I(U, B[i]).(gens(I(U)))) || error("ideals do not coincide on overlap")
          end
        end
      end
    end
    return I
  end
end

@doc raw"""
    PrimeIdealSheafFromChart

Type for sheaves of prime ideals ``P`` on a covered scheme ``X``
constructed from a prime ideal of the coordinate ring of a chart.
Essentially this is a scheme theoretic point.

For ``U`` an affine chart of ``X``, the ideal ``P(U)`` is computed using the gluings. 
The implementation is lazy.

# Examples
```jldoctest
julia> IP2 = projective_space(GF(7), [:x, :y, :z])
Projective space of dimension 2
  over prime field of characteristic 7
with homogeneous coordinates [x, y, z]

julia> X = covered_scheme(IP2)
Scheme
  over prime field of characteristic 7
with default covering
  described by patches
    1: affine 2-space
    2: affine 2-space
    3: affine 2-space
  in the coordinate(s)
    1: [(y//x), (z//x)]
    2: [(x//y), (z//y)]
    3: [(x//z), (y//z)]

julia> U, V, W = affine_charts(X);

julia> glue = default_covering(X)[U, V];

julia> y, z = gens(OO(U))
2-element Vector{FqMPolyRingElem}:
 (y//x)
 (z//x)

julia> P = ideal(OO(U), [y - z])
Ideal generated by
  (y//x) + 6*(z//x)

julia> PP = Oscar.PrimeIdealSheafFromChart(X, U, P)
Prime ideal sheaf on Scheme over GF(7) covered with 3 patches extended from Ideal ((y//x) + 6*(z//x)) on Affine 2-space

julia> PP(W)
Ideal generated by
  (y//z) + 6

```
"""
@attributes mutable struct PrimeIdealSheafFromChart{SpaceType, OpenType, OutputType,
                                                    RestrictionType
                                                   } <: AbsIdealSheaf{
                                                                      SpaceType, OpenType,
                                                                      OutputType, RestrictionType
                                                                     }
  X::AbsCoveredScheme
  U::AbsAffineScheme
  P::Ideal
  F::PreSheafOnScheme
  cheap_sub_ideals::WeakKeyIdDict{<:AbsAffineScheme, <:Ideal}

  function PrimeIdealSheafFromChart(
      X::AbsCoveredScheme,
      U::AbsAffineScheme,
      P::Ideal
    )
    @assert base_ring(P) === OO(U)
    @assert has_ancestor(x->any(y->y===x, affine_charts(X)), U) "the given affine scheme can not be matched with the affine charts of the covered scheme"

    Ipre = PreSheafOnScheme(X, 
                      OpenType=AbsAffineScheme, OutputType=Ideal,
                      RestrictionType=Map,
                      is_open_func=_is_open_func_for_schemes_without_affine_scheme_open_subscheme(X)
                     )
    cheap_sub_ideals = WeakKeyIdDict{AbsAffineScheme, Ideal}()
    I = new{typeof(X), AbsAffineScheme, Ideal, Map}(X, U, P, Ipre, cheap_sub_ideals)
    return I
  end
end

@doc raw"""
    SumIdealSheaf
    
Sum of two or more ideal sheaves.

# Examples
```jldoctest
julia> IP2 = projective_space(GF(7), [:x, :y, :z])
Projective space of dimension 2
  over prime field of characteristic 7
with homogeneous coordinates [x, y, z]

julia> X = covered_scheme(IP2)
Scheme
  over prime field of characteristic 7
with default covering
  described by patches
    1: affine 2-space
    2: affine 2-space
    3: affine 2-space
  in the coordinate(s)
    1: [(y//x), (z//x)]
    2: [(x//y), (z//y)]
    3: [(x//z), (y//z)]

julia> U, V, W = affine_charts(X);

julia> y, z = gens(OO(U));

julia> P1 = ideal(OO(U), [y - z]);

julia> PP1 = Oscar.PrimeIdealSheafFromChart(X, U, P1);

julia> P2 = ideal(OO(V), [OO(V)[1]]);

julia> PP2 = Oscar.PrimeIdealSheafFromChart(X, V, P2)
Prime ideal sheaf on Scheme over GF(7) covered with 3 patches extended from Ideal ((x//y)) on Affine 2-space

julia> II = PP1 + PP2
Sum of
  Prime ideal sheaf on scheme over GF(7) covered with 3 patches extended from ideal ((y//x) + 6*(z//x)) on affine 2-space
  Prime ideal sheaf on scheme over GF(7) covered with 3 patches extended from ideal ((x//y)) on affine 2-space

julia> typeof(II)
Oscar.SumIdealSheaf{CoveredScheme{FqField}, AbsAffineScheme, Ideal, Map}

julia> II(W)
Ideal generated by
  (y//z) + 6
  (x//z)

```
"""
@attributes mutable struct SumIdealSheaf{SpaceType, OpenType, OutputType,
                                         RestrictionType
                                        } <: AbsIdealSheaf{
                                                           SpaceType, OpenType,
                                                           OutputType, RestrictionType
                                                          }
  summands::Vector{<:AbsIdealSheaf}
  underlying_presheaf::AbsPreSheaf

  function SumIdealSheaf(
      summands::Vector{<:AbsIdealSheaf}
    )
    @assert !isempty(summands) "list of summands must not be empty"
    X = scheme(first(summands))
    @assert all(x->scheme(x) === X, summands)

    Ipre = PreSheafOnScheme(X,
                      OpenType=AbsAffineScheme, OutputType=Ideal,
                      RestrictionType=Map,
                      is_open_func=_is_open_func_for_schemes_without_affine_scheme_open_subscheme(X)
                     )
    I = new{typeof(X), AbsAffineScheme, Ideal, Map}(summands, Ipre)
    return I
  end
end

@doc raw"""
    ProductIdealSheaf

Product of two or more ideal sheaves.

# Examples
```jldoctest
julia> IP2 = projective_space(GF(7), [:x, :y, :z])
Projective space of dimension 2
  over prime field of characteristic 7
with homogeneous coordinates [x, y, z]

julia> X = covered_scheme(IP2)
Scheme
  over prime field of characteristic 7
with default covering
  described by patches
    1: affine 2-space
    2: affine 2-space
    3: affine 2-space
  in the coordinate(s)
    1: [(y//x), (z//x)]
    2: [(x//y), (z//y)]
    3: [(x//z), (y//z)]

julia> U, V, W = affine_charts(X);

julia> y, z = gens(OO(U));

julia> P1 = ideal(OO(U), [y - z]);

julia> PP1 = Oscar.PrimeIdealSheafFromChart(X, U, P1);

julia> P2 = ideal(OO(V), [OO(V)[1]]);

julia> PP2 = Oscar.PrimeIdealSheafFromChart(X, V, P2)
Prime ideal sheaf on Scheme over GF(7) covered with 3 patches extended from Ideal ((x//y)) on Affine 2-space

julia> II = PP1 * PP2
Product of
  Prime ideal sheaf on scheme over GF(7) covered with 3 patches extended from ideal ((y//x) + 6*(z//x)) on affine 2-space
  Prime ideal sheaf on scheme over GF(7) covered with 3 patches extended from ideal ((x//y)) on affine 2-space

julia> typeof(II)
Oscar.ProductIdealSheaf{CoveredScheme{FqField}, AbsAffineScheme, Ideal, Map}

julia> II(W)
Ideal generated by
  (x//z)*(y//z) + 6*(x//z)

```
"""
@attributes mutable struct ProductIdealSheaf{SpaceType, OpenType, OutputType,
                                             RestrictionType
                                            } <: AbsIdealSheaf{
                                                               SpaceType, OpenType,
                                                               OutputType, RestrictionType
                                                              }
  factors::Vector{<:AbsIdealSheaf}
  underlying_presheaf::AbsPreSheaf

  function ProductIdealSheaf(
      factors::Vector{<:AbsIdealSheaf}
    )
    @assert !isempty(factors) "list of summands must not be empty"
    X = scheme(first(factors))
    @assert all(x->scheme(x) === X, factors)

    Ipre = PreSheafOnScheme(X,
                      OpenType=AbsAffineScheme, OutputType=Ideal,
                      RestrictionType=Map,
                      is_open_func=_is_open_func_for_schemes_without_affine_scheme_open_subscheme(X)
                     )
    I = new{typeof(X), AbsAffineScheme, Ideal, Map}(factors, Ipre)
    return I
  end
end

@doc raw"""
    SimplifiedIdealSheaf
    
For a given `AbsIdealSheaf` `II` on an `AbsCoveredScheme` `X` this uses a 
heuristic to replace the generating set of `II(U)` by a hopefully smaller 
one on every affine chart of `X`.

# Examples
```jldoctes
julia> IP2 = projective_space(GF(7), [:x, :y, :z])
Projective space of dimension 2
  over prime field of characteristic 7
with homogeneous coordinates [x, y, z]

julia> X = covered_scheme(IP2)
Scheme
  over prime field of characteristic 7
with default covering
  described by patches
    1: affine 2-space
    2: affine 2-space
    3: affine 2-space
  in the coordinate(s)
    1: [(y//x), (z//x)]
    2: [(x//y), (z//y)]
    3: [(x//z), (y//z)]

julia> U, V, W = affine_charts(X);

julia> y, z = gens(OO(U));

julia> P1 = ideal(OO(U), [y - z]);

julia> PP1 = Oscar.PrimeIdealSheafFromChart(X, U, P1);

julia> P2 = ideal(OO(V), [OO(V)[1], OO(V)[2]]);

julia> PP2 = Oscar.PrimeIdealSheafFromChart(X, V, P2)
Prime ideal sheaf on Scheme over GF(7) covered with 3 patches extended from Ideal ((x//y), (z//y)) on Affine 2-space

julia> II = PP1 + PP2
Sum of 
  Prime ideal sheaf on scheme over GF(7) covered with 3 patches extended from ideal ((y//x) + 
  6*(z//x)) on affine 2-space
  Prime ideal sheaf on scheme over GF(7) covered with 3 patches extended from ideal ((x//y), (
  z//y)) on affine 2-space

julia> JJ = simplify(II)
Sheaf of ideals
  on scheme over GF(7) covered with 3 patches
    1: [(y//x), (z//x)]   affine 2-space
    2: [(x//y), (z//y)]   affine 2-space
    3: [(x//z), (y//z)]   affine 2-space
with restrictions
  1: Ideal (1)
  2: Ideal (1)
  3: Ideal (1)

julia> typeof(JJ)
Oscar.SimplifiedIdealSheaf{CoveredScheme{FqField}, AbsAffineScheme, Ideal, Map}

julia> II(W)
Ideal generated by
  (y//z) + 6
  1

julia> JJ(W)
Ideal generated by
  1

```
"""
@attributes mutable struct SimplifiedIdealSheaf{SpaceType, OpenType, OutputType,
                                             RestrictionType
                                            } <: AbsIdealSheaf{
                                                               SpaceType, OpenType,
                                                               OutputType, RestrictionType
                                                              }
  orig::AbsIdealSheaf
  underlying_presheaf::AbsPreSheaf

  function SimplifiedIdealSheaf(
      orig::AbsIdealSheaf
    )
    X = scheme(orig)

    Ipre = PreSheafOnScheme(X,
                      OpenType=AbsAffineScheme, OutputType=Ideal,
                      RestrictionType=Map,
                      is_open_func=_is_open_func_for_schemes_without_affine_scheme_open_subscheme(X)
                     )
    I = new{typeof(X), AbsAffineScheme, Ideal, Map}(orig, Ipre)
    return I
  end
end

@doc raw"""
    PullbackIdealSheaf
    
Given an morphism `f : X -> Y` of `AbsCoveredScheme`s and an ideal sheaf 
`II` on `Y`, this computes the pullback `f^* II` on `X`.

# Examples
```jldoctest
julia> IP2 = projective_space(GF(7), [:x, :y, :z])
Projective space of dimension 2
  over prime field of characteristic 7
with homogeneous coordinates [x, y, z]

julia> Y = covered_scheme(IP2)
Scheme
  over prime field of characteristic 7
with default covering
  described by patches
    1: affine 2-space
    2: affine 2-space
    3: affine 2-space
  in the coordinate(s)
    1: [(y//x), (z//x)]
    2: [(x//y), (z//y)]
    3: [(x//z), (y//z)]

julia> U, _ = affine_charts(Y);

julia> y, z = gens(OO(U));

julia> P1 = ideal(OO(U), [y, z]);

julia> PP1 = Oscar.PrimeIdealSheafFromChart(Y, U, P1);

julia> bl = blow_up(PP1);

julia> X = domain(bl);

julia> JJ = pullback(bl, PP1)
Sheaf of ideals
  on scheme over GF(7) covered with 4 patches
    1: [(s1//s0), (y//x)]   scheme(0)
    2: [(s0//s1), (z//x)]   scheme(0)
    3: [(x//y), (z//y)]     affine 2-space
    4: [(x//z), (y//z)]     affine 2-space
with restrictions
  1: Ideal ((y//x), (s1//s0)*(y//x))
  2: Ideal ((s0//s1)*(z//x), (z//x))
  3: Ideal (1)
  4: Ideal (1)

julia> typeof(JJ)
Oscar.PullbackIdealSheaf{CoveredScheme{FqField}, AbsAffineScheme, Ideal, Map}

```
"""
@attributes mutable struct PullbackIdealSheaf{SpaceType, OpenType, OutputType,
                                              RestrictionType
                                             } <: AbsIdealSheaf{
                                                                SpaceType, OpenType,
                                                                OutputType, RestrictionType
                                                               }
  f::AbsCoveredSchemeMorphism
  orig::AbsIdealSheaf
  Ipre::PreSheafOnScheme

  function PullbackIdealSheaf(
      f::AbsCoveredSchemeMorphism,
      orig::AbsIdealSheaf
    )
    X = domain(f)
    Y = codomain(f)
    @assert Y === scheme(orig)

    Ipre = PreSheafOnScheme(X,
                      OpenType=AbsAffineScheme, OutputType=Ideal,
                      RestrictionType=Map,
                      is_open_func=_is_open_func_for_schemes_without_affine_scheme_open_subscheme(X)
                     )
    I = new{typeof(X), AbsAffineScheme, Ideal, Map}(f, orig, Ipre)
    return I
  end
end

@doc raw"""
    RadicalOfIdealSheaf
    
Given an `AbsIdealSheaf` `II` on an `AbsCoveredScheme` `X`, this computes 
the sheaf associated to the radicals of the ideals on the charts.

# Examples
```jldoctest
julia> IP2 = projective_space(GF(7), [:x, :y, :z])
Projective space of dimension 2
  over prime field of characteristic 7
with homogeneous coordinates [x, y, z]

julia> Y = covered_scheme(IP2)
Scheme
  over prime field of characteristic 7
with default covering
  described by patches
    1: affine 2-space
    2: affine 2-space
    3: affine 2-space
  in the coordinate(s)
    1: [(y//x), (z//x)]
    2: [(x//y), (z//y)]
    3: [(x//z), (y//z)]

julia> U, _ = affine_charts(Y);

julia> y, z = gens(OO(U));

julia> P1 = ideal(OO(U), [y, z]);

julia> PP1 = Oscar.PrimeIdealSheafFromChart(Y, U, P1);

julia> PP1_squared = PP1^2;

julia> JJ = radical(PP1_squared);

julia> typeof(JJ)
Oscar.RadicalOfIdealSheaf{CoveredScheme{FqField}, AbsAffineScheme, Ideal, Map}

julia> PP1_squared(U)
Ideal generated by
  (y//x)^2
  (y//x)*(z//x)
  (z//x)^2

julia> JJ(U)
Ideal generated by
  (y//x)
  (y//x)*(z//x)
  (z//x)

```
"""
@attributes mutable struct RadicalOfIdealSheaf{SpaceType, OpenType, OutputType,
                                               RestrictionType
                                              } <: AbsIdealSheaf{
                                                                 SpaceType, OpenType,
                                                                 OutputType, RestrictionType
                                                                }
  orig::AbsIdealSheaf
  Ipre::PreSheafOnScheme

  function RadicalOfIdealSheaf(
      orig::AbsIdealSheaf
    )
    X = scheme(orig)
    Ipre = PreSheafOnScheme(X,
                      OpenType=AbsAffineScheme, OutputType=Ideal,
                      RestrictionType=Map,
                      is_open_func=_is_open_func_for_schemes_without_affine_scheme_open_subscheme(X)
                     )
    I = new{typeof(X), AbsAffineScheme, Ideal, Map}(orig, Ipre)
    set_attribute!(I, :is_radical=>true) # Set the attribute to prevent unnecessary computations
    return I
  end
end

@doc raw"""
    ToricIdealSheafFromCoxRingIdeal
    
This is the ideal sheaf associated to a `NormalToricVariety` `X` 
with `cox_ring` `S` and an ideal `I` of `S`.

# Examples
```jldoctest
julia> IP2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> S = cox_ring(IP2)
Multivariate polynomial ring in 3 variables over QQ graded by
  x1 -> [1]
  x2 -> [1]
  x3 -> [1]

julia> x, y, z = gens(S);

julia> I = ideal(S, x^2 - y*z);

julia> II = ideal_sheaf(IP2, I)
Sheaf of ideals
  on normal toric variety
with restrictions
  1: Ideal (x_1_1^2 - x_2_1)
  2: Ideal (x_1_2*x_2_2 - 1)
  3: Ideal (x_1_3^2 - x_2_3)

julia> typeof(II)
Oscar.ToricIdealSheafFromCoxRingIdeal{NormalToricVariety, AbsAffineScheme, Ideal, Map}

```
"""
@attributes mutable struct ToricIdealSheafFromCoxRingIdeal{SpaceType, OpenType, OutputType,
                                                           RestrictionType
                                                          } <: AbsIdealSheaf{
                                                                             SpaceType, OpenType,
                                                                             OutputType, RestrictionType
                                                                            }

  X::NormalToricVariety
  I::MPolyIdeal
  underlying_presheaf::PreSheafOnScheme

  function ToricIdealSheafFromCoxRingIdeal(X::NormalToricVariety, I::MPolyIdeal)
    @assert cox_ring(X) === base_ring(I) "incompatible rings"
    Ipre = PreSheafOnScheme(X,
                      OpenType=AbsAffineScheme, OutputType=Ideal,
                      RestrictionType=Map,
                      is_open_func=_is_open_func_for_schemes_without_affine_scheme_open_subscheme(X)
                     )
    I = new{typeof(X), AbsAffineScheme, Ideal, Map}(X, I, Ipre)
    return I
  end
end

########################################################################
# Singular locus ideal sheaf
########################################################################

@doc raw"""
    SingularLocusIdealSheaf
    
This is the (radical) ideal sheaf for the singular locus of an 
`AbsCoveredScheme` `X`.

# Examples
```jldoctest
julia> IP2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> S = cox_ring(IP2)
Multivariate polynomial ring in 3 variables over QQ graded by
  x1 -> [1]
  x2 -> [1]
  x3 -> [1]

julia> x, y, z = gens(S);

julia> I = ideal(S, x^3 - y*z^2);

julia> II = ideal_sheaf(IP2, I);

julia> X, inc_X = sub(II);

julia> JJ = Oscar.ideal_sheaf_of_singular_locus(X)
Sheaf of ideals
  on scheme over QQ covered with 3 patches
    1: [x_1_1, x_2_1]   scheme(x_1_1^3 - x_2_1)
    2: [x_1_2, x_2_2]   scheme(x_1_2^2*x_2_2 - 1)
    3: [x_1_3, x_2_3]   scheme(x_1_3^3 - x_2_3^2)
with restrictions
  1: Ideal (1)
  2: Ideal (1)
  3: Ideal (x_1_3, x_2_3)

julia> typeof(JJ)
Oscar.SingularLocusIdealSheaf{CoveredScheme{QQField}, AbsAffineScheme, Ideal, Map}

```
"""
@attributes mutable struct SingularLocusIdealSheaf{SpaceType, OpenType, OutputType,
                                         RestrictionType
                                        } <: AbsIdealSheaf{
                                                           SpaceType, OpenType,
                                                           OutputType, RestrictionType
                                                          }
  focus::AbsIdealSheaf
  non_radical_ideals::IdDict{<:AbsAffineScheme, <:Ideal}
  underlying_presheaf::AbsPreSheaf

  function SingularLocusIdealSheaf(
      X::AbsCoveredScheme;
      focus::AbsIdealSheaf=zero_ideal_sheaf(X)
    )

    Ipre = PreSheafOnScheme(X,
                      OpenType=AbsAffineScheme, OutputType=Ideal,
                      RestrictionType=Map,
                      is_open_func=_is_open_func_for_schemes_without_affine_scheme_open_subscheme(X)
                     )
    I = new{typeof(X), AbsAffineScheme, Ideal, Map}(focus, IdDict{AbsAffineScheme, Ideal}(), Ipre)
    return I
  end
end

underlying_presheaf(I::SingularLocusIdealSheaf) = I.underlying_presheaf
focus(I::SingularLocusIdealSheaf) = I.focus


