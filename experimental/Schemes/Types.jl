
export EmptyScheme
export AbsSpec, Spec
export ClosedEmbedding
export AbsGlueing, Glueing
export SimpleGlueing
export AbsProjectiveScheme, ProjectiveScheme
export ProjectiveSchemeMor
export Covering, CoveringMorphism
export AbsCoveredScheme, CoveredScheme
export AbsCoveredSchemeMorphism, CoveredSchemeMorphism
export VarietyFunctionField, VarietyFunctionFieldElem
export IdealSheaf




########################################################################
# Special Type for closed embeddings of affine schemes                 #
########################################################################
@Markdown.doc """
    ClosedEmbedding{DomainType, CodomainType, PullbackType}

A closed embedding ``f : X → Y`` of affine schemes ``X = Spec(S)``
into ``Y = Spec(R)`` such that ``S ≅ R/I`` via ``f`` for some
ideal ``I ⊂ R``.
"""
@attributes mutable struct ClosedEmbedding{DomainType,
                                           CodomainType,
                                           PullbackType
                                          }<:AbsSpecMor{DomainType,
                                                        CodomainType,
                                                        PullbackType,
                                                        ClosedEmbedding,
                                                        Nothing
                                                       }
  inc::SpecMor{DomainType, CodomainType, PullbackType}
  I::Ideal
  U::SpecOpen

  function ClosedEmbedding(X::AbsSpec, I::Ideal)
    base_ring(I) == OO(X) || error("ideal does not belong to the correct ring")
    Y = subscheme(X, I)
    inc = SpecMor(Y, X, hom(OO(X), OO(Y), gens(OO(Y))))
    return new{typeof(Y), typeof(X), pullback_type(inc)}(inc, I)
  end
  function ClosedEmbedding(f::SpecMor, I::Ideal; check::Bool=true)
    Y = domain(f)
    X = codomain(f)
    base_ring(I) == OO(X) || error("ideal does not belong to the correct ring")
    if check
      Y == subscheme(X, I)
      pullback(f).(gens(OO(X))) == gens(OO(Y))
    end
    return new{typeof(Y), typeof(X), pullback_type(f)}(f, I)
  end
end

########################################################################
# Abstract glueings for affine schemes                                 #
########################################################################
abstract type AbsGlueing{LeftSpecType<:AbsSpec,
                         RightSpecType<:AbsSpec,
                         LeftOpenType<:Scheme,
                         RightOpenType<:Scheme,
                         LeftMorType<:Hecke.Map,
                         RightMorType<:Hecke.Map
                        } end

########################################################################
# Concrete type for general glueings                                   #
########################################################################
@Markdown.doc """
    Glueing{SpecType<:Spec, OpenType<:SpecOpen, MorType<:SpecOpenMor}

Glueing of two affine schemes ``X ↩ U ≅ V ↪ Y`` along open subsets
``U ⊂ X`` and ``V ⊂ Y via some isomorphism ``φ : U → V``.
"""
@attributes mutable struct Glueing{
                                   LeftSpecType<:AbsSpec,
                                   RightSpecType<:AbsSpec,
                                   LeftOpenType<:SpecOpen,
                                   RightOpenType<:SpecOpen,
                                   LeftMorType<:SpecOpenMor,
                                   RightMorType<:SpecOpenMor
                                  } <: AbsGlueing{
                                   LeftSpecType,
                                   RightSpecType,
                                   LeftOpenType,
                                   RightOpenType,
                                   LeftMorType,
                                   RightMorType
                                  }
  X::LeftSpecType
  Y::RightSpecType
  U::LeftOpenType
  V::RightOpenType
  f::LeftMorType # f : U → V
  g::RightMorType

  function Glueing(
      X::AbsSpec, Y::AbsSpec, f::SpecOpenMor, g::SpecOpenMor; check::Bool=true
    )
    ambient_scheme(domain(f)) === X || error("the domain of the glueing morphism is not an open subset of the first argument")
    ambient_scheme(codomain(f)) === Y || error("the codomain of the glueing morphism is not an open subset of the second argument")
    if check
      (domain(f) === codomain(g) &&
      domain(g) ===  codomain(f)) || error("maps can not be isomorphisms")
      compose(f, g) == identity_map(domain(f)) || error("glueing maps are not inverse of each other")
      compose(g, f) == identity_map(domain(g)) || error("glueing maps are not inverse of each other")
    end
    return new{typeof(X), typeof(Y),
               typeof(domain(f)), typeof(domain(g)),
               typeof(f), typeof(g)
              }(X, Y, domain(f), domain(g), f, g)
  end
end

########################################################################
# Special type for simple glueings of affine schemes along principal
# open subsets
#
# SimpleGlueing is for glueings X ↩ U ≅ V ↪ Y along principal
# open subsets U ⊂ X and V ⊂ Y along identifications f : U ↔ V : g.
# For general glueings it can not be guaranteed to have this setup,
# but it is a situation often encountered and with significant
# simplification of underlying algorithms in the background.
# Hence, the special type.
########################################################################
@attributes mutable struct SimpleGlueing{LST<:AbsSpec,
                                         RST<:AbsSpec,
                                         LOT<:PrincipalOpenSubset,
                                         ROT<:PrincipalOpenSubset,
                                         LMT<:AbsSpecMor,
                                         RMT<:AbsSpecMor
                                        } <: AbsGlueing{LST, RST, LOT, ROT, LMT, RMT}
  X::LST
  Y::RST
  U::LOT
  V::ROT
  f::LMT
  g::RMT

  function SimpleGlueing(
      X::AbsSpec, Y::AbsSpec,
      f::AbsSpecMor{<:PrincipalOpenSubset},
      g::AbsSpecMor{<:PrincipalOpenSubset};
      check::Bool=true
    )
    U = domain(f)
    V = domain(g)
    X === ambient_scheme(U) && Y === ambient_scheme(V) || error("schemes are not compatible")
    domain(f) === codomain(g) && domain(g) === codomain(f) || error("maps are not compatible")
    if check
      is_identity_map(compose(f, g)) || error("maps are not inverse to each other")
      is_identity_map(compose(g, f)) || error("maps are not inverse to each other")
    end
    set_attribute!(f, :inverse, g)
    set_attribute!(g, :inverse, f)
    return new{typeof(X), typeof(Y),
               typeof(U), typeof(V),
               typeof(f), typeof(g)
              }(X, Y, U, V, f, g)
  end
end

########################################################################
# Abstract projective schemes                                          #
########################################################################
abstract type AbsProjectiveScheme{BaseRingType, RingType} <: Scheme{BaseRingType} end

########################################################################
# Concrete type for projective schemes                                 #
########################################################################
@Markdown.doc """
    ProjectiveScheme{CoeffRingType, CoeffRingElemType, RingType, RingElemType}

Closed subschemes ``X ⊂ ℙʳ(A)`` of projective space of `fiber_dimension` ``r``
over a ring of coefficients ``A`` of type `CoeffRingType` with elements of
type `CoeffRingElemType`. The subscheme ``X`` is given by means of a homogeneous
ideal ``I`` in the graded ring ``A[s₀,…,sᵣ]`` and the latter is of type
`RingType` with elements of type `RingElemType`.
"""
@attributes mutable struct ProjectiveScheme{CoeffRingType, CoeffRingElemType, RingType, RingElemType} <: AbsProjectiveScheme{CoeffRingType, RingType}
  A::CoeffRingType	# the base ring
  r::Int	# the relative dimension
  S::RingType   # A[s₀,…,sᵣ]
  I::MPolyIdeal{RingElemType} # generators for the defining ideal

  # fields used for caching
  C::Scheme # The affine cone of this scheme.
  Y::Scheme # the base scheme
  projection_to_base::SchemeMor
  homog_coord::Vector # the homogeneous coordinates as functions on the affine cone

  function ProjectiveScheme(S::MPolyRing_dec)
    all(x->(total_degree(x) == 1), gens(S)) || error("ring is not standard graded")
    n = ngens(S)-1
    A = coefficient_ring(S)
    I = ideal(S, [zero(S)])
    return new{typeof(A), elem_type(A), typeof(S), elem_type(S)}(A, n, S, I)
  end

  function ProjectiveScheme(S::MPolyRing_dec, I::MPolyIdeal{T}) where {T<:RingElem}
    for f in gens(I)
      parent(f) == S || error("elements do not belong to the correct ring")
    end
    all(x->(total_degree(x) == 1), gens(S)) || error("ring is not standard graded")
    n = ngens(S)-1
    A = coefficient_ring(S)
    return new{typeof(A), elem_type(A), typeof(S), elem_type(S)}(A, n, S, I)
  end

  function ProjectiveScheme(Q::MPolyQuo{MPolyElem_dec{T, AbstractAlgebra.Generic.MPoly{T}}}) where {T}
    all(x->(total_degree(x) == 1), gens(S)) || error("ring is not standard graded")
    S = base_ring(Q)
    A = coefficient_ring(S)
    I = gens(modulus(Q))
    r = ngens(S)-1
    return new{typeof(A), elem_type(A), typeof(S), elem_type(S)}(A, r, S, I)
  end
end

########################################################################
# Morphisms of projective schemes                                      #
########################################################################
@Markdown.doc """
    ProjectiveSchemeMor

A morphism of projective schemes

    ℙˢ(B)     ℙʳ(A)
      ∪         ∪
      P    →    Q
      ↓         ↓
   Spec(B) → Spec(A)

given by means of a commutative diagram of homomorphisms of
graded rings

  A[v₀,…,vᵣ] → B[u₀,…,uₛ]
      ↑            ↑
      A      →     B

If no morphism `A → B` of the base rings is specified, then
both ``P`` and ``Q`` are assumed to be defined in relative projective
space over the same ring with the identity on the base.
"""
mutable struct ProjectiveSchemeMor{
    DomainType<:ProjectiveScheme,
    CodomainType<:ProjectiveScheme,
    PullbackType<:Hecke.Map,
    BaseMorType
  } <: SchemeMor{DomainType, CodomainType,
                 ProjectiveSchemeMor,
                 BaseMorType
                }
  domain::DomainType
  codomain::CodomainType
  pullback::PullbackType

  #fields for caching
  map_on_base_schemes::SchemeMor
  map_on_affine_cones::SchemeMor

  ### Simple morphism of projective schemes over the same base scheme
  function ProjectiveSchemeMor(
      P::DomainType,
      Q::CodomainType,
      f::PullbackType;
      check::Bool=true
    ) where {DomainType<:ProjectiveScheme, CodomainType<:ProjectiveScheme, PullbackType<:Map}
    T = ambient_ring(P)
    S = ambient_ring(Q)
    (S === domain(f) && T === codomain(f)) || error("pullback map incompatible")
    if check
      #TODO: Check map on ideals (not available yet)
    end
    return new{DomainType, CodomainType, PullbackType, Nothing}(P, Q, f)
  end

  ### complicated morphisms over a non-trivial morphism of base schemes
  function ProjectiveSchemeMor(
      P::DomainType,
      Q::CodomainType,
      f::PullbackType,
      h::BaseMorType;
      check::Bool=true
    ) where {DomainType<:ProjectiveScheme,
             CodomainType<:ProjectiveScheme,
             PullbackType<:Map,
             BaseMorType<:SchemeMor
            }
    T = ambient_ring(P)
    S = ambient_ring(Q)
    (S === domain(f) && T === codomain(f)) || error("pullback map incompatible")
    pbh = pullback(h)
    OO(domain(h)) == coefficient_ring(T) || error("base scheme map not compatible")
    OO(codomain(h)) == coefficient_ring(S) || error("base scheme map not compatible")
    if check
      T(pbh(one(OO(codomain(h))))) == f(S(one(OO(codomain(h))))) == one(T) || error("maps not compatible")
      coefficient_map(f) == pbh || error("maps not compatible")
    end
    return new{DomainType, CodomainType, PullbackType, BaseMorType}(P, Q, f, h)
  end
end

########################################################################
# Coverings for covered schemes                                        #
########################################################################
@Markdown.doc """
    Covering

A covering of a scheme ``X`` by affine patches ``Uᵢ`` which are glued
along isomorphisms ``gᵢⱼ : Uᵢ⊃ Vᵢⱼ →  Vⱼᵢ ⊂ Uⱼ``.

**Note:** The distinction between the different affine patches of the scheme
is made from their hashes. Thus, an affine scheme must not appear more than once
in any covering!
"""
mutable struct Covering{BaseRingType}
  patches::Vector{<:AbsSpec} # the basic affine patches of X
  glueings::IdDict{Tuple{<:AbsSpec, <:AbsSpec}, <:AbsGlueing} # the glueings of the basic affine patches
  affine_refinements::IdDict{<:AbsSpec, <:Vector{<:Tuple{<:SpecOpen, Vector{<:RingElem}}}} # optional lists of refinements
      # of the basic affine patches.
      # These are stored as pairs (U, a) where U is a 'trivial' SpecOpen,
      # meaning that its list of hypersurface equation (f₁,…,fᵣ) has empty
      # intersection in the basic affine patch X and hence satisfies
      # some equality 1 ≡ a₁⋅f₁ + a₂⋅f₂ + … + aᵣ⋅fᵣ on X.
      # Since the coefficients aᵢ of this equality are crucial for computations,
      # we store them in an extra tuple.

  # fields for caching
  glueing_graph::Graph{Undirected}
  transition_graph::Graph{Undirected}
  edge_dict::Dict{Tuple{Int, Int}, Int}

  function Covering(
      patches::Vector{<:AbsSpec},
      glueings::IdDict{Tuple{<:AbsSpec, <:AbsSpec}, <:AbsGlueing};
      check::Bool=true,
      affine_refinements::IdDict{
          <:AbsSpec,
          <:Vector{<:Tuple{<:SpecOpen, <:Vector{<:RingElem}}}
         }=IdDict{AbsSpec, Vector{Tuple{SpecOpen, Vector{RingElem}}}}()
    )
    n = length(patches)
    n > 0 || error("can not glue the empty scheme")
    kk = coefficient_ring(ambient_ring(patches[1]))
    for i in 2:n
      kk == coefficient_ring(base_ring(OO(patches[i]))) || error("schemes are not defined over the same base ring")
    end
    # Check that no patch appears twice
    for i in 1:n-1
      for j in i+1:n
        patches[i] === patches[j] && error("affine schemes must not appear twice among the patches")
      end
    end
    for (X, Y) in keys(glueings)
      X in patches || error("glueings are not compatible with the patches")
      Y in patches || error("glueings are not compatible with the patches")
      if haskey(glueings, (Y, X))
        if check
          inverse(glueings[(X, Y)]) == glueings[(Y, X)] || error("glueings are not inverse of each other")
        end
      else
        glueings[(Y, X)] = inverse(glueings[(X, Y)])
      end
    end

    # check the affine refinements
    for U in keys(affine_refinements)
      for (V, a) in affine_refinements[U]
        ambient_scheme(V) == U && error("the ambient scheme of the refinement of X must be X")
        U in patches && error("the ambient scheme of the refinement can not be found in the affine patches")
        if check
          isone(OO(U)(sum([c*g for (c, g) in zip(a, gens(U))]))) || error("the patch $V does not cover $U")
        end
      end
    end
    return new{base_ring_type(patches[1])}(patches, glueings, affine_refinements)
  end

  ### the empty covering
  function Covering(kk::Ring)
    return new{typeof(kk)}(Vector{AbsSpec}(), IdDict{Tuple{AbsSpec, AbsSpec}, AbsGlueing}(),
                           IdDict{AbsSpec, Vector{Tuple{SpecOpen, Vector{RingElem}}}}())
  end
end

########################################################################
# Morphisms of coverings                                               #
########################################################################
@Markdown.doc """
    CoveringMorphism{SpecType<:Spec, CoveringType<:Covering, SpecMorType<:SpecMor}

A morphism ``f : C → D`` of two coverings. For every patch ``U`` of ``C`` this
provides a map `f[U']` of type `SpecMorType` from ``U' ⊂ U`` to
some patch `codomain(f[U])` in `D` for some affine patches ``U'`` covering ``U``.

**Note:** For two affine patches ``U₁, U₂ ⊂ U`` the codomains of `f[U₁]` and `f[U₂]`
do not need to coincide! However, given the glueings in `C` and `D`, all affine maps
have to coincide on their overlaps.
"""
mutable struct CoveringMorphism{DomainType<:Covering, CodomainType<:Covering, BaseMorType}
  domain::DomainType
  codomain::CodomainType
  morphisms::IdDict{<:AbsSpec, <:AbsSpecMor} # on a patch X of the domain covering, this
                                         # returns the morphism φ : X → Y to the corresponding
                                         # patch Y of the codomain covering.

  function CoveringMorphism(
      dom::DomainType,
      cod::CodomainType,
      mor::IdDict{<:AbsSpec, <:AbsSpecMor};
      check::Bool=true
    ) where {
             DomainType<:Covering,
             CodomainType<:Covering
            }
    # TODO: check domain/codomain compatibility
    # TODO: if check is true, check that all morphisms glue and that the domain patches
    # cover the basic patches of `dom`.
    for U in keys(mor)
      U in dom || error("patch $U of the map not found in domain")
      codomain(mor[U]) in cod || error("codomain patch not found")
    end
    # check that the whole domain is covered
    for U in basic_patches(dom)
      if !haskey(mor, U)
        !haskey(affine_refinements(dom), U) || error("patch $U of the domain not covered")
        found = false
        for (V, a) in affine_refinements(dom)[U]
          all(x->(haskey(mor, x)), affine_patches(V)) && (found = true)
        end
        !found && error("patch $U of the domain not covered")
      end
    end
    return new{DomainType, CodomainType, Nothing}(dom, cod, mor)
  end
end

########################################################################
# Abstract type for covered schemes                                    #
########################################################################
abstract type AbsCoveredScheme{BaseRingType} <: Scheme{BaseRingType} end

########################################################################
# A minimal implementation of AbsCoveredScheme                         #
########################################################################
@Markdown.doc """
    mutable struct CoveredScheme{
      CoveringType<:Covering,
      CoveringMorphismType<:CoveringMorphism
    }

A covered scheme ``X`` given by means of at least one covering
of type `CoveringType`.

A scheme may posess several coverings which are partially ordered
by refinement. Such refinements are special instances of `CoveringMorphism`

    ρ : C1 → C2

where for each patch ``U`` in `C1` the inclusion map ``ρ[U] : U → V``
into the corresponding patch ``V`` of `C2` is an open embedding for which
both ``𝒪(U)`` and ``𝒪(V)`` have the same `base_ring` (so that they can be
canonically compared).
"""
@attributes mutable struct CoveredScheme{BaseRingType} <: AbsCoveredScheme{BaseRingType}
  coverings::Vector{<:Covering}
  refinements::Dict{<:Tuple{<:Covering, <:Covering}, <:CoveringMorphism}
  refinement_graph::Graph{Directed}
  kk::BaseRingType

  default_covering::Covering

  function CoveredScheme(coverings::Vector{<:Covering},
      refinements::Dict{Tuple{<:Covering, <:Covering}, <:CoveringMorphism}
    )
    # TODO: Check whether the refinements form a connected graph.
    BaseRingType = base_ring_type(coverings[1])
    all(x->(base_ring_type(x) == BaseRingType), coverings) || error("coverings are not compatible")
    X = new{BaseRingType}(coverings, refinements)
    X.default_covering = X.coverings[1]
    X.kk = base_ring(patches(coverings[1])[1])
    return X
  end
  function CoveredScheme(kk::Ring)
    res = new{typeof(kk)}()
    res.kk = kk
    return res
  end
end

########################################################################
# Morphisms of covered schemes                                         #
########################################################################
abstract type AbsCoveredSchemeMorphism{
    DomainType<:CoveredScheme,
    CodomainType<:CoveredScheme,
    BaseMorphismType,
    CoveredSchemeMorphismType
   } <: SchemeMor{DomainType, CodomainType, CoveredSchemeMorphismType, BaseMorphismType}
end

########################################################################
# Concrete minimal type for morphisms of covered schemes               #
########################################################################
@attributes mutable struct CoveredSchemeMorphism{
    DomainType<:AbsCoveredScheme,
    CodomainType<:AbsCoveredScheme,
    BaseMorphismType
   } <: AbsCoveredSchemeMorphism{
                                 DomainType,
                                 CodomainType,
                                 CoveredSchemeMorphism,
                                 BaseMorphismType
                                }
  X::DomainType
  Y::CodomainType
  f::CoveringMorphism

  function CoveredSchemeMorphism(
      X::DomainType,
      Y::CodomainType,
      f::CoveringMorphism{<:Any, <:Any, BaseMorType};
      check::Bool=true
    ) where {
             DomainType<:CoveredScheme,
             CodomainType<:CoveredScheme,
             BaseMorType
            }
    domain(f) in coverings(X) || error("covering not found in domain")
    codomain(f) in coverings(Y) || error("covering not found in codomain")
    return new{DomainType, CodomainType, BaseMorType}(X, Y, f)
  end
end

########################################################################
# Rational functions on irreducible varieties                          #
########################################################################

mutable struct VarietyFunctionField{BaseRingType<:Field,
                                    FracFieldType<:AbstractAlgebra.Generic.FracField,
                                    CoveredSchemeType<:AbsCoveredScheme,
                                    SpecType<:AbsSpec
                                   } <: Field
  kk::BaseRingType
  X::CoveredSchemeType
  U::SpecType  # representative patch to represent rational functions
  KK::FracFieldType

  function VarietyFunctionField(
      X::AbsCoveredScheme;
      check::Bool=true,
      representative_patch::AbsSpec=default_covering(X)[1]
    )
    check && (is_irreducible(X) || error("variety is not irreducible"))
    representative_patch in default_covering(X) || error("representative patch not found")
    KK = FractionField(ambient_ring(representative_patch))
    kk = base_ring(X)
    return new{typeof(kk), typeof(KK), typeof(X), typeof(representative_patch)}(kk, X, representative_patch, KK)
  end
end

########################################################################
# Elements of VarietyFunctionFields                                    #
########################################################################
mutable struct VarietyFunctionFieldElem{FracType<:AbstractAlgebra.Generic.Frac,
                                        ParentType<:VarietyFunctionField
                                       }
  KK::ParentType
  f::FracType

  function VarietyFunctionFieldElem(
      KK::VarietyFunctionField,
      f::AbstractAlgebra.Generic.Frac;
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
@Markdown.doc """
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

@Markdown.doc """
    PreSheafOnScheme

A basic minimal implementation of the interface for `AbsPreSheaf`; to be used internally.
"""
@attributes mutable struct PreSheafOnScheme{SpaceType, OpenType, OutputType, RestrictionType,
                                       IsOpenFuncType, ProductionFuncType,
                                       RestrictionFuncType
                                      } <: AbsPreSheaf{
                                       SpaceType, OpenType,
                                       OutputType, RestrictionType
                                      }
  X::SpaceType

  # caches
  obj_cache::IdDict{<:OpenType, <:OutputType} # To cache values that have already been computed
  res_cache::IdDict{<:Tuple{<:OpenType, <:OpenType}, <:RestrictionType} # To cache already computed restrictions

  # production functions for new objects
  is_open_func::IsOpenFuncType # To check whether one set is open in the other
  production_func::ProductionFuncType # To produce ℱ(U) for U ⊂ X
  restriction_func::RestrictionFuncType

  function PreSheafOnScheme(X::Scheme, production_func::Any, restriction_func::Any;
      OpenType=AbsSpec, OutputType=Any, RestrictionType=Any,
      is_open_func::Any=is_open_embedding
    )
    return new{typeof(X), OpenType, OutputType, RestrictionType,
               typeof(is_open_func), typeof(production_func), typeof(restriction_func)
              }(X, IdDict{OpenType, OutputType}(),
                IdDict{Tuple{OpenType, OpenType}, RestrictionType}(),
                is_open_func, production_func, restriction_func
               )
  end
end

########################################################################
# The structure sheaf of affine and covered schemes                    #
########################################################################
@Markdown.doc """
    StructureSheafOfRings <: AbsPreSheaf

On an `AbsCoveredScheme` ``X`` this returns the sheaf ``𝒪`` of rings of
regular functions on ``X``.

Note that due to technical reasons, the admissible open subsets are restricted
to the following:
 * `U::AbsSpec` among the `basic_patches` of the `default_covering` of `X`;
 * `U::PrincipalOpenSubset` with `ambient_scheme(U)` in the `basic_patches` of the `default_covering` of `X`;
 * `W::SpecOpen` with `ambient_scheme(W)` in the `basic_patches` of the `default_covering` of `X`.

One can call the restriction maps of ``𝒪`` across charts, implicitly using the
identifications given by the glueings in the `default_covering`.
"""
@attributes mutable struct StructureSheafOfRings{SpaceType, OpenType, OutputType,
                                          RestrictionType, ProductionFuncType,
                                          RestrictionFuncType,
                                          PreSheafType
                                         } <: AbsPreSheaf{
                                          SpaceType, OpenType,
                                          OutputType, RestrictionType
                                         }
  OO::PreSheafType

  ### Structure sheaf on affine schemes
  function StructureSheafOfRings(X::AbsSpec)
    function is_open_func(U::AbsSpec, V::AbsSpec)
      return is_subset(V, X) && is_open_embedding(U, V) # Note the restriction to subsets of X
    end
    function production_func(U::AbsSpec)
      return OO(U)
    end
    function restriction_func(V::AbsSpec, OV::Ring, U::AbsSpec, OU::Ring)
      return hom(OV, OU, gens(OU), check=false) # check=false assures quicker computation
    end

    R = PreSheafOnScheme(X, production_func, restriction_func,
                    OpenType=AbsSpec, OutputType=Ring,
                    RestrictionType=Hecke.Map,
                    is_open_func=is_open_func
                   )
    return new{typeof(X), Union{AbsSpec, SpecOpen}, Ring, Hecke.Map,
               typeof(production_func), typeof(restriction_func),
               typeof(R)}(R)
  end

  ### Structure sheaf on covered schemes
  function StructureSheafOfRings(X::AbsCoveredScheme)

    ### Checks for open containment.
    #
    # We allow the following cases:
    #
    #  * U::PrincipalOpenSubset in W===ambient_scheme(U) in the basic charts of X
    #  * U::PrincipalOpenSubset ⊂ V::PrincipalOpenSubset with ambient_scheme(U) === ambient_scheme(V) in the basic charts of X
    #  * U::PrincipalOpenSubset ⊂ V::PrincipalOpenSubset with ambient_scheme(U) != ambient_scheme(V) both in the basic charts of X
    #    and U and V contained in the glueing domains of their ambient schemes
    #  * U::AbsSpec ⊂ U::AbsSpec in the basic charts of X
    #  * U::AbsSpec ⊂ X for U in the basic charts
    #  * U::PrincipalOpenSubset ⊂ X with ambient_scheme(U) in the basic charts of X
    #  * W::SpecOpen ⊂ X with ambient_scheme(U) in the basic charts of X
    function is_open_func(U::PrincipalOpenSubset, V::PrincipalOpenSubset)
      C = default_covering(X)
      A = ambient_scheme(U)
      A in C || return false
      B = ambient_scheme(V)
      B in C || return false
      if A === B
        is_subset(U, V) || return false
      else
        G = C[A, B] # Get the glueing
        f, g = glueing_morphisms(G)
        is_subset(U, domain(f)) || return false
        is_subset(V, domain(g)) || return false
        gU = preimage(g, U)
        is_subset(gU, V) || return false
      end
      return true
    end
    function is_open_func(U::PrincipalOpenSubset, Y::AbsCoveredScheme)
      return Y === X && ambient_scheme(U) in default_covering(X)
    end
    function is_open_func(U::AbsSpec, Y::AbsCoveredScheme)
      return Y === X && U in default_covering(X)
    end
    function is_open_func(Z::AbsCoveredScheme, Y::AbsCoveredScheme)
      return X === Y === Z
    end
    function is_open_func(U::AbsSpec, V::AbsSpec)
      U in default_covering(X) || return false
      V in default_covering(X) || return false
      G = default_covering(X)[U, V]
      return issubset(U, glueing_domains(G)[1])
    end
    function is_open_func(U::PrincipalOpenSubset, V::AbsSpec)
      V in default_covering(X) || return false
      ambient_scheme(U) === V && return true
      W = ambient_scheme(U)
      W in default_covering(X) || return false
      G = default_covering(X)[W, V]
      return is_subset(U, glueing_domains(G)[1])
    end
#    function is_open_func(U::PrincipalOpenSubset, V::PrincipalOpenSubset)
#      ambient_scheme(V) in default_covering(X) || return false
#      ambient_scheme(U) === ambient_scheme(V) && return issubset(U, V)
#      W = ambient_scheme(U)
#      W in default_covering(X) || return false
#      G = default_covering(X)[W, ambient_scheme(V)]
#      preV = preimage(glueing_morphisms(G)[1], V)
#      return is_subset(U, preV)
#    end
    function is_open_func(W::SpecOpen, Y::AbsCoveredScheme)
      return Y === X && ambient_scheme(W) in default_covering(X)
    end
    function is_open_func(W::SpecOpen, V::AbsSpec)
      V in default_covering(X) || return false
      ambient_scheme(W) === V && return true
      U = ambient_scheme(W)
      U in default_covering(X) || return false
      G = default_covering(X)[U, V]
      return is_subset(W, glueing_domains(G)[1])
    end
    function is_open_func(W::SpecOpen, V::PrincipalOpenSubset)
      PW = ambient_scheme(W)
      PV = ambient_scheme(V)
      PW in default_covering(X) || return false
      PV in default_covering(X) || return false
      if PW === PV
        return issubset(W, V)
        #return all(x->(issubset(x, V)), affine_patches(W))
      else
        G = default_covering(X)[PW, PV]
        preV = preimage(glueing_morphisms(G)[1], V)
        return issubset(W, preV)
      end
    end
    function is_open_func(W::SpecOpen, V::SpecOpen)
      PW = ambient_scheme(W)
      PV = ambient_scheme(V)
      PW in default_covering(X) || return false
      PV in default_covering(X) || return false
      if PW === PV
        return issubset(W, V)
        #return all(x->(issubset(x, V)), affine_patches(W))
      else
        G = default_covering(X)[PW, PV]
        preV = preimage(glueing_morphisms(G)[1], V)
        return issubset(W, preV)
      end
    end
    function is_open_func(U::AbsSpec, W::SpecOpen)
      U in default_covering(X) || return false
      if U === ambient_scheme(W)
        # in this case W must be equal to U
        return issubset(W, U)
        #return one(OO(U)) in complement_ideal(W)
      else
        G = default_covering(X)[ambient_scheme(W), U]
        issubset(U, glueing_domains(G)[2]) || return false
        preU = preimage(glueing_morphisms(G)[1], U)
        return issubset(preU, W)
      end
    end
    function is_open_func(U::PrincipalOpenSubset, W::SpecOpen)
      ambient_scheme(U) in default_covering(X) || return false
      if ambient_scheme(U) === ambient_scheme(W)
        # in this case W must be equal to U
        return issubset(W, U)
        #return one(OO(U)) in complement_ideal(W)
      else
        G = default_covering(X)[ambient_scheme(W), ambient_scheme(U)]
        issubset(U, glueing_domains(G)[2]) || return false
        preU = preimage(glueing_morphisms(G)[1], U)
        return issubset(preU, W)
      end
    end

    ### Production of the rings of regular functions; to be cached
    function production_func(U::AbsSpec)
      return OO(U)
    end
    function production_func(U::SpecOpen)
      return OO(U)
    end

    ### Production of the restriction maps; to be cached
    function restriction_func(V::AbsSpec, OV::Ring, U::AbsSpec, OU::Ring)
      V === U || error("basic affine patches must be the same")
      return identity_map(OV)
    end
    function restriction_func(V::AbsSpec, OV::Ring, U::PrincipalOpenSubset, OU::Ring)
      if ambient_scheme(U) === V
        return hom(OV, OU, gens(OU), check=false)
      else
        W = ambient_scheme(U)
        G = default_covering(X)[V, W]
        f, g = glueing_morphisms(G)
        pbg = pullback(g)
        function rho_func(x::RingElem)
          parent(x) == OV || error("element does not belong to the correct domain")
          return restrict(pbg(domain(pbg)(x)), U) # should probably be tuned to avoid checks.
        end
        return hom(OV, OU, rho_func.(gens(OV)), check=false)
      end
      error("arguments are not valid")
    end
    function restriction_func(V::PrincipalOpenSubset, OV::Ring, U::AbsSpec, OU::Ring)
      if ambient_scheme(V) === U
        function rho_func(a::RingElem)
          parent(a) === OV || error("element does not belong to the correct ring")
          # We may assume that all denominators admissible in V are
          # already units in OO(U)
          return OU(lifted_numerator(a))*inv(OU(lifted_denominator(a)))
        end
        return hom(OV, OU, rho_func.(gens(OV)), check=false)
      else
        G = default_covering(X)[ambient_scheme(V), U]
        W1, W2 = glueing_domains(G)
        f, g = glueing_morphisms(G)
        function rho_func2(a::RingElem)
          parent(a) === OV || error("element does not belong to the correct ring")
          return restrict(pullback(g)(OO(W1)(a)), U)
        end
        return hom(OV, OU, rho_func2.(gens(OV)), check=false)
      end
    end
    function restriction_func(V::PrincipalOpenSubset, OV::Ring, U::PrincipalOpenSubset, OU::Ring)
      A = ambient_scheme(V)
      if A === ambient_scheme(U)
        return hom(OV, OU, gens(OU), check=false)
      else
        B = ambient_scheme(U)
        G = default_covering(X)[A, B]
        f, g = glueing_morphisms(G)
        function rho_func(x::RingElem)
          parent(x) == OV || error("input not valid")
          y = pullback(g)(OO(codomain(g))(x))
          return restrict(pullback(g)(OO(codomain(g))(x)), U)
        end
        return hom(OV, OU, rho_func.(gens(OV)), check=false)
      end
      error("arguments are invalid")
    end
    function restriction_func(V::AbsSpec, OV::Ring, W::SpecOpen, OW::Ring)
      V in default_covering(X) || return false
      ambient_scheme(W) in default_covering(X) || return false
      if V === ambient_scheme(W)
        return MapFromFunc(x->(OW(x)), OV, OW)
      else
        G = default_covering(X)[V, ambient_scheme(W)]
        f, g = glueing_morphisms(G)
        function rho_func(a::RingElem)
          parent(a) === OV || error("element does not belong to the correct ring")
          return restrict(pullback(g)(OO(domain(f))(a)), W)
        end
        return MapFromFunc(rho_func, OV, OW)
      end
    end
    function restriction_func(V::PrincipalOpenSubset, OV::Ring, W::SpecOpen, OW::Ring)
      if ambient_scheme(V) === ambient_scheme(W)
        function rho_func(a::RingElem)
          parent(a) === OV || error("element does not belong to the correct ring")
          return OW(a)
        end
        return MapFromFunc(rho_func, OV, OW)
      else
        G = default_covering(X)(ambient_scheme(V), ambient_scheme(W))
        f, g = glueing_morphisms(G)
        VG = intersect(V, domain(f))
        preV = preimage(g, VG)
        gres = restriction(g, preV, VG, check=false)
        inc = inclusion_morphism(W, preV)
        function rho_func2(a::RingElem)
          parent(a) === OV || error("element does not belong to the correct ring")
          return pullback(inc)(pullback(gres)(OO(preV)(a)))
        end
        return MapFromFunc(rho_func2, OV, OW)
      end
    end
    function restriction_func(V::SpecOpen, OV::Ring, W::SpecOpen, OW::Ring)
      if ambient_scheme(V) === ambient_scheme(W)
        inc = inclusion_morphism(W, V)
        return MapFromFunc(pullback(inc), OV, OW)
      else
        G = default_covering(X)[ambient_scheme(V), ambient_scheme(W)]
        f, g = glueing_morphisms(G)
        VG = intersect(V, domain(f))
        inc0 = inclusion_morphism(VG, V)
        preV = preimage(g, VG)
        gres = restrict(g, preV, VG, check=false)
        inc = inclusion_morphism(W, preV)
        return MapFromFunc(x->(pullback(inc)(pullback(gres)(pullback(inc0)(x)))),
                           OV, OW)
      end
    end

    R = PreSheafOnScheme(X, production_func, restriction_func,
                      OpenType=Union{AbsSpec, SpecOpen}, OutputType=Ring,
                      RestrictionType=Hecke.Map,
                      is_open_func=is_open_func
                     )
    return new{typeof(X), Union{AbsSpec, SpecOpen}, Ring, Hecke.Map,
               typeof(production_func), typeof(restriction_func),
               typeof(R)}(R)
  end
end

########################################################################
# Ideal sheaves on covered schemes                                     #
########################################################################
@Markdown.doc """
    IdealSheaf <: AbsPreSheaf

A sheaf of ideals ``ℐ`` on an `AbsCoveredScheme` ``X``.

Note that due to technical reasons, the admissible open subsets are restricted
to the following:
 * `U::AbsSpec` among the `basic_patches` of the `default_covering` of `X`;
 * `U::PrincipalOpenSubset` with `ambient_scheme(U)` in the `basic_patches` of the `default_covering` of `X`.

One can call the restriction maps of ``ℐ`` across charts, implicitly using the
identifications given by the glueings in the `default_covering`.
"""
@attributes mutable struct IdealSheaf{SpaceType, OpenType, OutputType,
                                      RestrictionType, ProductionFuncType,
                                      RestrictionFuncType,
                                      PreSheafType
                                     } <: AbsPreSheaf{
                                                      SpaceType, OpenType,
                                                      OutputType, RestrictionType
                                                     }
  ID::IdDict{AbsSpec, Ideal} # the ideals on the basic patches of the default covering
  OOX::StructureSheafOfRings # the structure sheaf on X
  I::PreSheafType # the underlying presheaf of ideals for caching

  ### Ideal sheaves on covered schemes
  function IdealSheaf(X::AbsCoveredScheme, ID::IdDict{AbsSpec, Ideal};
      check::Bool=true
    )
    OOX = StructureSheafOfRings(X)

    ### Checks for open containment.
    #
    # We allow the following cases:
    #
    #  * U::PrincipalOpenSubset in W===ambient_scheme(U) in the basic charts of X
    #  * U::PrincipalOpenSubset ⊂ V::PrincipalOpenSubset with ambient_scheme(U) === ambient_scheme(V) in the basic charts of X
    #  * U::PrincipalOpenSubset ⊂ V::PrincipalOpenSubset with ambient_scheme(U) != ambient_scheme(V) both in the basic charts of X
    #    and U and V contained in the glueing domains of their ambient schemes
    #  * U::AbsSpec ⊂ U::AbsSpec in the basic charts of X
    #  * U::AbsSpec ⊂ X for U in the basic charts
    #  * U::PrincipalOpenSubset ⊂ X with ambient_scheme(U) in the basic charts of X
    function is_open_func(U::PrincipalOpenSubset, V::PrincipalOpenSubset)
      C = default_covering(X)
      A = ambient_scheme(U)
      A in C || return false
      B = ambient_scheme(V)
      B in C || return false
      if A === B
        is_subset(U, V) || return false
      else
        G = C[A, B] # Get the glueing
        f, g = glueing_morphisms(G)
        is_subset(U, domain(f)) || return false
        is_subset(V, domain(g)) || return false
        gU = preimage(g, U)
        is_subset(gU, V) || return false
      end
      return true
    end
    function is_open_func(U::PrincipalOpenSubset, Y::AbsCoveredScheme)
      return Y === X && ambient_scheme(U) in default_covering(X)
    end
    function is_open_func(U::AbsSpec, Y::AbsCoveredScheme)
      return Y === X && U in default_covering(X)
    end
    function is_open_func(Z::AbsCoveredScheme, Y::AbsCoveredScheme)
      return X === Y === Z
    end
    function is_open_func(U::AbsSpec, V::AbsSpec)
      U in default_covering(X) || return false
      V in default_covering(X) || return false
      G = default_covering(X)[U, V]
      return issubset(U, glueing_domains(G)[1])
    end
    function is_open_func(U::PrincipalOpenSubset, V::AbsSpec)
      V in default_covering(X) || return false
      ambient_scheme(U) === V && return true
      W = ambient_scheme(U)
      W in default_covering(X) || return false
      G = default_covering(X)[W, V]
      return is_subset(U, glueing_domains(G)[1])
    end

    ### Production of the rings of regular functions; to be cached
    function production_func(U::AbsSpec)
      haskey(ID, U) && return ID[U]
      # The ideal sheaf has to be provided on at one dense
      # open subset of every connected component.
      # Otherwise, the ideal sheaf is given by the unit
      # ideals.
      for G in glueings(default_covering(X))
        A, B = patches(G)
        Asub, Bsub = glueing_domains(G)
        if A === U && haskey(ID, B) && is_dense(Asub)
          Z = subscheme(B, ID[B])
          f, _ = glueing_morphisms(G)
          pZ = preimage(f, Z)
          ZU = closure(pZ, U)
          ID[U] = ideal(OO(U), gens(saturated_ideal(modulus(OO(ZU)))))
          return ID[U]
        end
      end
      return ideal(OO(U), one(OO(U)))
    end
    function production_func(U::PrincipalOpenSubset)
      V = ambient_scheme(U)
      IV = production_func(V)
      rho = OOX(V, U)
      IU = ideal(OO(U), rho.(gens(IV)))
      return IU
    end

    ### Production of the restriction maps; to be cached
    function restriction_func(V::AbsSpec, IV::Ideal, U::AbsSpec, IU::Ideal)
      return OOX(V, U) # This does not check containment of the arguments
                       # in the ideal. But this is not a parent check and
                       # hence expensive, so we might want to not do that.
    end

    Ipre = PreSheafOnScheme(X, production_func, restriction_func,
                      OpenType=AbsSpec, OutputType=Ideal,
                      RestrictionType=Hecke.Map,
                      is_open_func=is_open_func
                     )
    I = new{typeof(X), AbsSpec, Ideal, Hecke.Map,
               typeof(production_func), typeof(restriction_func),
               typeof(Ipre)}(ID, OOX, Ipre)
    if check
      # Check that all ideal sheaves are compatible on the overlaps.
      # TODO: eventually replace by a check that on every basic
      # affine patch, the ideal sheaf can be inferred from what is
      # given on one dense open subset.
      C = default_covering(X)
      for U in basic_patches(default_covering(X))
        for V in basic_patches(default_covering(X))
          G = C[U, V]
          A, B = glueing_domains(G)
          for i in 1:ngens(A)
            I(A[i]) == ideal(OOX(A[i]), I(V, A[i]).(gens(I(V)))) || error("ideals do not coincide on overlap")
          end
          for i in 1:ngens(B)
            I(B[i]) == ideal(OOX(B[i]), I(U, B[i]).(gens(I(U)))) || error("ideals do not coincide on overlap")
          end
        end
      end
    end
    return I
  end
end
