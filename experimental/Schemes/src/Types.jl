export AbsPreSheaf
export AbsProjectiveScheme
export EmptyScheme
export IdealSheaf
export ProjectiveScheme
export ProjectiveSchemeMor
export VarietyFunctionField
export VarietyFunctionFieldElem


########################################################################
# Rational functions on irreducible varieties                          #
########################################################################

mutable struct VarietyFunctionField{BaseRingType<:Field,
                                    FracFieldType<:AbstractAlgebra.Generic.FracField,
                                    CoveredSchemeType<:AbsCoveredScheme,
                                    AffineSchemeType<:AbsAffineScheme
                                   } <: Field
  kk::BaseRingType
  X::CoveredSchemeType
  U::AffineSchemeType  # representative patch to represent rational functions
  KK::FracFieldType

  function VarietyFunctionField(
      X::AbsCoveredScheme;
      check::Bool=true,
      representative_patch::AbsAffineScheme=default_covering(X)[1]
    )
    @check is_irreducible(X) "variety is not irreducible"
    representative_patch in default_covering(X) || error("representative patch not found")
    KK = fraction_field(ambient_coordinate_ring(representative_patch))
    kk = base_ring(X)
    return new{typeof(kk), typeof(KK), typeof(X), typeof(representative_patch)}(kk, X, representative_patch, KK)
  end
end

########################################################################
# Elements of VarietyFunctionFields                                    #
########################################################################
mutable struct VarietyFunctionFieldElem{FracType<:AbstractAlgebra.Generic.FracFieldElem,
                                        ParentType<:VarietyFunctionField
                                       }
  KK::ParentType
  f::FracType

  function VarietyFunctionFieldElem(
      KK::VarietyFunctionField,
      f::AbstractAlgebra.Generic.FracFieldElem;
      check::Bool=true
    )
    representative_field(KK) == parent(f) || error("element does not have the correct parent")
    return new{typeof(f), typeof(KK)}(KK, f)
  end

  function VarietyFunctionFieldElem(
      KK::VarietyFunctionField,
      a::RingElem, b::RingElem;
      check::Bool=true
    )
    R = parent(a)
    R == parent(b) || error("parent rings not compatible")
    R == base_ring(representative_field(KK))
    f = representative_field(KK)(a, b)
    return new{typeof(f), typeof(KK)}(KK, f)
  end
end

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

  # production functions for new objects
  is_open_func::Function # To check whether one set is open in the other
  production_func::Union{Function, Nothing} # To produce ℱ(U) for U ⊂ X
  restriction_func::Union{Function, Nothing} # To produce the restriction maps ℱ(U) → ℱ(V) for V ⊂ U ⊂ X open

  function PreSheafOnScheme(X::Scheme, 
      production_func::Union{Function, Nothing}=nothing, 
      restriction_func::Union{Function, Nothing}=nothing;
      OpenType=AbsAffineScheme, OutputType=Any, RestrictionType=Any,
      is_open_func::Any=is_open_embedding
    )
    return new{typeof(X), OpenType, OutputType, RestrictionType,
              }(X, IdDict{OpenType, OutputType}(),
                #IdDict{Tuple{OpenType, OpenType}, RestrictionType}(),
                is_open_func, production_func, restriction_func
               )
  end
end

########################################################################
# Simplified Spectra                                                   #
########################################################################
@attributes mutable struct SimplifiedAffineScheme{BaseRingType, RingType<:Ring} <: AbsAffineScheme{BaseRingType, RingType}
  X::AbsAffineScheme
  Y::AbsAffineScheme
  f::AbsAffineSchemeMor
  g::AbsAffineSchemeMor

  function SimplifiedAffineScheme(X::AbsAffineScheme, Y::AbsAffineScheme, f::AbsAffineSchemeMor, g::AbsAffineSchemeMor;
      check::Bool=true
    )
    domain(f) === X || error("map is not compatible")
    codomain(f) === Y || error("map is not compatible")
    domain(g) === Y || error("map is not compatible")
    codomain(g) === X || error("map is not compatible")

    @check is_identity_map(compose(f, g)) && is_identity_map(compose(g, f)) "maps are not inverse to each other"

    result = new{typeof(base_ring(X)), typeof(OO(X))}(X, Y)
    # We need to rewrap the identification maps so that the (co-)domains match
    fwrap = morphism(result, Y, pullback(f), check=false)
    gwrap = morphism(Y, result, pullback(g), check=false)
    set_attribute!(fwrap, :inverse, gwrap)
    set_attribute!(gwrap, :inverse, fwrap)
    result.f = fwrap
    result.g = gwrap
    return result 
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

    function production_func(F::AbsPreSheaf, U::AbsAffineScheme)
      return OO(U)
    end
    function restriction_func(F::AbsPreSheaf, V::AbsAffineScheme, U::AbsAffineScheme)
      OU = F(U) # assumed to be previously cached
      OV = F(V) # same as above
      return hom(OV, OU, gens(OU), check=false) # check=false assures quicker computation
    end

    R = PreSheafOnScheme(X, production_func, restriction_func,
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

########################################################################
# Closed embeddings                                                    #
########################################################################

@attributes mutable struct CoveredClosedEmbedding{
    DomainType<:AbsCoveredScheme,
    CodomainType<:AbsCoveredScheme,
    BaseMorphismType
   } <: AbsCoveredSchemeMorphism{
                                 DomainType,
                                 CodomainType,
                                 BaseMorphismType,
                                 CoveredSchemeMorphism
                                }
  f::CoveredSchemeMorphism
  I::AbsIdealSheaf

  function CoveredClosedEmbedding(
      X::DomainType,
      Y::CodomainType,
      f::CoveringMorphism{<:Any, <:Any, MorphismType, BaseMorType};
      check::Bool=true,
      ideal_sheaf::AbsIdealSheaf=IdealSheaf(Y, f, check=check)
    ) where {
             DomainType<:AbsCoveredScheme,
             CodomainType<:AbsCoveredScheme,
             MorphismType<:ClosedEmbedding,
             BaseMorType
            }
    ff = CoveredSchemeMorphism(X, Y, f; check)
    if has_decomposition_info(codomain(f))
      for U in patches(domain(f))
        floc = f[U]
        phi = pullback(floc)
        V = codomain(floc)
        g = Vector{elem_type(OO(V))}(decomposition_info(codomain(f))[V])
        set_decomposition_info!(domain(f), U, Vector{elem_type(OO(U))}(phi.(g)))
      end
    end
    #all(x->(x isa ClosedEmbedding), values(morphisms(f))) || error("the morphisms on affine patches must be `ClosedEmbedding`s")
    return new{DomainType, CodomainType, BaseMorType}(ff, ideal_sheaf)
  end
end

########################################################################
# PushforwardIdealSheaf                                                #
########################################################################
@attributes mutable struct PushforwardIdealSheaf{SpaceType, OpenType, OutputType,
                                                 RestrictionType
                                                } <: AbsIdealSheaf{
                                                                   SpaceType, OpenType,
                                                                   OutputType, RestrictionType
                                                                  }
  f::CoveredClosedEmbedding
  orig::AbsIdealSheaf
  Ipre::PreSheafOnScheme

  function PushforwardIdealSheaf(
      f::CoveredClosedEmbedding,
      orig::AbsIdealSheaf
    )
    X = domain(f)
    Y = codomain(f)
    @assert X === scheme(orig)

    Ipre = PreSheafOnScheme(Y,
                      OpenType=AbsAffineScheme, OutputType=Ideal,
                      RestrictionType=Map,
                      is_open_func=_is_open_func_for_schemes_without_affine_scheme_open_subscheme(X)
                     )
    I = new{typeof(Y), AbsAffineScheme, Ideal, Map}(f, orig, Ipre)
    return I
  end
end

########################################################################
# Morphisms from rational functions                                    #
########################################################################
@doc raw"""
    MorphismFromRationalFunctions{DomainType<:AbsCoveredScheme, CodomainType<:AbsCoveredScheme} 

A lazy type for a dominant morphism ``φ : X → Y`` of `AbsCoveredScheme`s which is given 
by a set of rational functions ``a₁,…,aₙ`` in the fraction field of the `base_ring`
of ``𝒪(U)`` for one of the dense open `affine_chart`s ``U`` of ``X``. 
The ``aᵢ`` represent the pullbacks of the coordinates (`gens`) of some 
`affine_chart` ``V`` of the codomain ``Y`` under this map. 
```jldoctest
julia> IP1 = covered_scheme(projective_space(QQ, [:s, :t]))
Scheme
  over rational field
with default covering
  described by patches
    1: affine 1-space
    2: affine 1-space
  in the coordinate(s)
    1: [(t//s)]
    2: [(s//t)]

julia> IP2 = projective_space(QQ, [:x, :y, :z]);

julia> S = homogeneous_coordinate_ring(IP2);

julia> x, y, z = gens(S);

julia> IPC, inc_IPC = sub(IP2, ideal(S, [x^2 - y*z]));

julia> C = covered_scheme(IPC);

julia> U = first(affine_charts(IP1))
Spectrum
  of multivariate polynomial ring in 1 variable (t//s)
    over rational field

julia> V = first(affine_charts(C))
Spectrum
  of quotient
    of multivariate polynomial ring in 2 variables (y//x), (z//x)
      over rational field
    by ideal (-(y//x)*(z//x) + 1)

julia> t = first(gens(OO(U)))
(t//s)

julia> Phi = MorphismFromRationalFunctions(IP1, C, U, V, [t//one(t), 1//t]);

julia> realizations = Oscar.realize_on_patch(Phi, U);

julia> realizations[3]
Affine scheme morphism
  from [(t//s)]          AA^1
  to   [(x//z), (y//z)]  scheme((x//z)^2 - (y//z))
given by the pullback function
  (x//z) -> (t//s)
  (y//z) -> (t//s)^2

```
"""
@attributes mutable struct MorphismFromRationalFunctions{DomainType<:AbsCoveredScheme, 
                                       CodomainType<:AbsCoveredScheme
                                      } <: AbsCoveredSchemeMorphism{DomainType, CodomainType, 
                                                                    MorphismFromRationalFunctions, Nothing}
  domain::DomainType
  codomain::CodomainType
  domain_covering::Covering
  codomain_covering::Covering
  domain_chart::AbsAffineScheme
  codomain_chart::AbsAffineScheme
  coord_imgs::Vector{<:FieldElem}
  run_internal_checks::Bool

  ### Various fields for caching
  patch_representatives::IdDict{<:AbsAffineScheme, <:Tuple{<:AbsAffineScheme, <:Vector{<:FieldElem}}}
  realizations::IdDict{<:AbsAffineScheme, <:Vector{<:AffineSchemeMor}}
  realization_previews::IdDict{<:Tuple{<:AbsAffineScheme, <:AbsAffineScheme}, <:Vector{<:FieldElem}}
  maximal_extensions::IdDict{<:Tuple{<:AbsAffineScheme, <:AbsAffineScheme}, <:Vector{<:AbsAffineSchemeMor}}
  cheap_realizations::IdDict{<:Tuple{<:AbsAffineScheme, <:AbsAffineScheme}, <:AbsAffineSchemeMor}
  full_realization::CoveredSchemeMorphism

  function MorphismFromRationalFunctions(
      X::AbsCoveredScheme, Y::AbsCoveredScheme, 
      U::AbsAffineScheme, V::AbsAffineScheme,
      a::Vector{<:FieldElem};
      check::Bool=true,
      domain_covering::Covering=default_covering(X),
      codomain_covering::Covering=default_covering(Y)
    )
    @check is_irreducible(X) "domain must be irreducible"
    @check is_irreducible(Y) "codomain must be irreducible"
    @check dim(Y) <= dim(X) "cannot be dominant"
    #_find_chart(U, default_covering(X)) !== nothing || error("patch not found in domain")
    #_find_chart(V, default_covering(Y)) !== nothing || error("patch not found in codomain")
    any(x->x===U, patches(default_covering(X))) || error("patch not found in domain")
    any(x->x===V, patches(default_covering(Y))) || error("patch not found in codomain")
    F = parent(first(a))
    R = base_ring(F)
    all(x->parent(x)===F, a) || error("coordinate images must be elements of the same field")
    R === ambient_coordinate_ring(U) || error("images of pullback of the coordinates do not live in the correct ring")
    patch_repr = IdDict{AbsAffineScheme, Tuple{AbsAffineScheme, Vector{FieldElem}}}()
    patch_repr[U] = (V, a)
    realizations = IdDict{AbsAffineScheme, Vector{AffineSchemeMor}}()
    realization_previews = IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, Vector{FieldElem}}()
    maximal_extensions = IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, Vector{AbsAffineSchemeMor}}()
    cheap_realizations = IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, AbsAffineSchemeMor}()
    return new{typeof(X), typeof(Y)}(X, Y, domain_covering, codomain_covering, 
                                     U, V, a, check, patch_repr, realizations, 
                                     realization_previews, maximal_extensions,
                                     cheap_realizations
                                    )
  end
end

