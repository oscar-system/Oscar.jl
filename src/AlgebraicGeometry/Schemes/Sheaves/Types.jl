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

 * `restriction_map(F, U, V)` for *admissible* open subsets ``V ⊂ U ⊂ X``: This returns the restriction map ``ρ : ℱ(U) → ℱ(V)``.
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

@doc """raw
    PrimeIdealSheafFromChart

Type for sheaves of prime ideals ``P`` on a covered scheme ``X``
constructed from a prime ideal of the coordinate ring of a chart.
Essentially this is a scheme theoretic point.

For ``U`` an affine chart of ``X``, the ideal ``P(U)`` is computed using the gluings. 
The implementation is lazy.
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
    set_attribute!(I, :is_radical=>true) # Some methods might be blind to is_radical and only want to check `is_known_to_be_radical` via attributes. Setting this makes sure they get it. 
    return I
  end
end

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


