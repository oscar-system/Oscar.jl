export EmptyScheme
export AbsProjectiveScheme, ProjectiveScheme
export ProjectiveSchemeMor
export VarietyFunctionField, VarietyFunctionFieldElem
export AbsPreSheaf
export IdealSheaf


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
@attributes mutable struct ProjectiveSchemeMor{
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
  base_ring_morphism::Hecke.Map

  #fields for caching
  map_on_base_schemes::SchemeMor
  map_on_affine_cones::SchemeMor

  ### Simple morphism of projective schemes over the same base scheme
  function ProjectiveSchemeMor(
      P::DomainType,
      Q::CodomainType,
      f::PullbackType;
      check::Bool=true
    ) where {DomainType<:ProjectiveScheme, 
             CodomainType<:ProjectiveScheme, 
             PullbackType<:Map
            }
    T = ambient_coordinate_ring(P)
    S = ambient_coordinate_ring(Q)
    (S === domain(f) && T === codomain(f)) || error("pullback map incompatible")
    if check
      #TODO: Check map on ideals (not available yet)
    end
    return new{DomainType, CodomainType, PullbackType, Nothing}(P, Q, f)
  end
  
  ### Morphisms with an underlying base change
  function ProjectiveSchemeMor(
      P::DomainType,
      Q::CodomainType,
      f::PullbackType;
      check::Bool=true
    ) where {DomainType<:ProjectiveScheme, 
             CodomainType<:ProjectiveScheme, 
             PullbackType<:MPolyAnyMap{<:Any, <:Any, <:Map}
            }
    T = ambient_coordinate_ring(P)
    S = ambient_coordinate_ring(Q)
    (S === domain(f) && T === codomain(f)) || error("pullback map incompatible")
    if check
      #TODO: Check map on ideals (not available yet)
    end
    return new{DomainType, CodomainType, PullbackType, Nothing}(P, Q, f, coefficient_map(f))
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
    T = ambient_coordinate_ring(P)
    S = ambient_coordinate_ring(Q)
    (S === domain(f) && T === codomain(f)) || error("pullback map incompatible")
    pbh = pullback(h)
    OO(domain(h)) == coefficient_ring(T) || error("base scheme map not compatible")
    OO(codomain(h)) == coefficient_ring(S) || error("base scheme map not compatible")
    if check
      T(pbh(one(OO(codomain(h))))) == f(S(one(OO(codomain(h))))) == one(T) || error("maps not compatible")
      coefficient_map(f) == pbh || error("maps not compatible")
    end
    return new{DomainType, CodomainType, PullbackType, BaseMorType}(P, Q, f, coefficient_map(f), h)
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
    KK = FractionField(ambient_coordinate_ring(representative_patch))
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
  production_func::Function # To produce ℱ(U) for U ⊂ X
  restriction_func::Function # To produce the restriction maps ℱ(U) → ℱ(V) for V ⊂ U ⊂ X open

  function PreSheafOnScheme(X::Scheme, production_func::Any, restriction_func::Any;
      OpenType=AbsSpec, OutputType=Any, RestrictionType=Any,
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
@attributes mutable struct SimplifiedSpec{BaseRingType, RingType} <: AbsSpec{BaseRingType, RingType} 
  X::AbsSpec
  Y::AbsSpec
  f::AbsSpecMor
  g::AbsSpecMor

  function SimplifiedSpec(X::AbsSpec, Y::AbsSpec, f::AbsSpecMor, g::AbsSpecMor;
      check::Bool=true
    )
    domain(f) === X || error("map is not compatible")
    codomain(f) === Y || error("map is not compatible")
    domain(g) === Y || error("map is not compatible")
    codomain(g) === X || error("map is not compatible")

    if check
      is_identity_map(compose(f, g)) && is_identity_map(compose(g, f)) || error("maps are not inverse to each other")
    end

    result = new{typeof(base_ring(X)), typeof(OO(X))}(X, Y)
    # We need to rewrap the identification maps so that the (co-)domains match
    fwrap = SpecMor(result, Y, pullback(f))
    gwrap = SpecMor(Y, result, pullback(g))
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
                                                 RestrictionType
                                                } <: AbsPreSheaf{
                                                                 SpaceType, OpenType,
                                                                 OutputType, RestrictionType
                                                                }
  OO::PreSheafOnScheme

  ### Structure sheaf on affine schemes
  function StructureSheafOfRings(X::AbsSpec)
    function is_open_func(U::AbsSpec, V::AbsSpec)
      return is_subset(V, X) && is_open_embedding(U, V) # Note the restriction to subsets of X
    end

    function production_func(F::AbsPreSheaf, U::AbsSpec)
      return OO(U)
    end
    function restriction_func(F::AbsPreSheaf, V::AbsSpec, U::AbsSpec)
      OU = F(U) # assumed to be previously cached
      OV = F(V) # same as above
      return hom(OV, OU, gens(OU), check=false) # check=false assures quicker computation
    end

    R = PreSheafOnScheme(X, production_func, restriction_func,
                    OpenType=AbsSpec, OutputType=Ring,
                    RestrictionType=Hecke.Map,
                    is_open_func=is_open_func
                   )
    return new{typeof(X), Union{AbsSpec, SpecOpen}, Ring, Hecke.Map}(R)
  end

  ### Structure sheaf on covered schemes
  function StructureSheafOfRings(X::AbsCoveredScheme)

    ### Production of the rings of regular functions; to be cached
    function production_func(F::AbsPreSheaf, U::AbsSpec)
      return OO(U)
    end
    function production_func(F::AbsPreSheaf, U::SpecOpen)
      return OO(U)
    end

    ### Production of the restriction maps; to be cached
    function restriction_func(F::AbsPreSheaf, V::AbsSpec, U::AbsSpec)
      V === U || error("basic affine patches must be the same")
      return identity_map(OO(V))
    end
    function restriction_func(
        F::AbsPreSheaf, 
        V::AbsSpec, 
        U::Union{<:PrincipalOpenSubset, <:SimplifiedSpec}
      )
      OV = F(V) # Assumed to be cached or produced on the fly.
      OU = F(U) # Same as above.
      incU = _flatten_open_subscheme(U, default_covering(X))
      #incU, dU = _find_chart(U, default_covering(X))
      U_flat = codomain(incU)
      W = ambient_scheme(U_flat)
      if W === V
        return pullback(compose(incU, inclusion_morphism(U_flat)))
      else
        G = default_covering(X)[V, W]
        f, g = glueing_morphisms(G)
        pbg = pullback(g)
        function rho_func(x::RingElem)
          parent(x) === OV || error("element does not belong to the correct domain")
          return pullback(incU)(restrict(pbg(domain(pbg)(x)), U_flat)) # should probably be tuned to avoid checks.
        end
        return hom(OV, OU, rho_func.(gens(OV)), check=false)
      end
      error("arguments are not valid")
    end

    function restriction_func(F::AbsPreSheaf, 
        V::Union{<:PrincipalOpenSubset, <:SimplifiedSpec},
        U::AbsSpec
      )
      OV = F(V)
      OU = F(U) 
      incV = _flatten_open_subscheme(V, default_covering(X))
      W = ambient_scheme(codomain(incV))
      V_direct = domain(incV)
      if W === U
        # By virtue of the checks in _is_open_func we must have V isomorphic to U.
        phi = pullback(inverse(incV))
        psi = hom(OO(V_direct), OU, gens(OU))
        return hom(OV, OU, psi.(phi.(gens(OV))))
        ### deprecated code below;
        # kept for the moment because of possible incompatibilities with glueings 
        # along SpecOpens.
        function rho_func(a::RingElem)
          parent(a) === OV || error("element does not belong to the correct ring")
          # We may assume that all denominators admissible in V are
          # already units in OO(U)
          return OU(lifted_numerator(a))*inv(OU(lifted_denominator(a)))
        end
        return hom(OV, OU, rho_func.(gens(OV)), check=false)
      else
        G = default_covering(X)[W, U]
        W1, W2 = glueing_domains(G)
        f, g = glueing_morphisms(G)
        g_res = restrict(g, U, V_direct)
        return pullback(compose(g_res, inverse(incV)))
        ### deprecated code below; see comment above
        function rho_func2(a::RingElem)
          parent(a) === OV || error("element does not belong to the correct ring")
          return restrict(pullback(g)(OO(W1)(a)), U)
        end
        return hom(OV, OU, rho_func2.(gens(OV)), check=false)
      end
    end
    function restriction_func(F::AbsPreSheaf, 
        V::Union{<:PrincipalOpenSubset, <:SimplifiedSpec},
        U::Union{<:PrincipalOpenSubset, <:SimplifiedSpec}
      )
      OV = F(V)
      OU = F(U)
      inc_U_flat = _flatten_open_subscheme(U, default_covering(X))
      inc_V_flat = _flatten_open_subscheme(V, default_covering(X))
      A = ambient_scheme(codomain(inc_U_flat))
      B = ambient_scheme(codomain(inc_V_flat))
      U_flat = codomain(inc_U_flat)
      V_flat = codomain(inc_V_flat)

      if A === B
        return hom(OV, OU, 
                   pullback(inc_U_flat).(pullback(inclusion_morphism(U_flat, V_flat)).(pullback(inverse(inc_V_flat)).(gens(OV)))), check=false
                  )
      else
        G = default_covering(X)[A, B]
        f, g = glueing_morphisms(G)
        VV_flat = intersect(V_flat, codomain(f))
        VU = preimage(f, VV_flat)
        fres = restrict(f, VU, VV_flat)
        inc_V_flat_inv = inverse(inc_V_flat)
        function rho_func(x::RingElem)
          parent(x) === OV || error("input not valid")
          y = pullback(inverse(inc_V_flat))(x)
          y = restrict(y, VV_flat)
          y = pullback(fres)(y)
          y = restrict(y, U_flat)
          return pullback(inc_U_flat)(y)
        end
        return hom(OV, OU, rho_func.(gens(OV)), check=false)
      end
      error("arguments are invalid")
    end

    function restriction_func(F::AbsPreSheaf, V::AbsSpec, W::SpecOpen)
      OV = F(V)
      OW = F(W)
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

    ### cleaned up until here ###
    # We do not make SpecOpen compatible with the tree structures, yet. 
    # All SpecOpen's are hence required to have an ambient_scheme on the top level. 

    function restriction_func(F::AbsPreSheaf, V::PrincipalOpenSubset, W::SpecOpen)
      error("method not implemented at the moment")
      OV = F(V)
      OW = F(W)
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
    function restriction_func(F::AbsPreSheaf, V::SpecOpen, W::SpecOpen)
      OV = F(V)
      OW = F(W)
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
                      is_open_func=_is_open_func_for_schemes(X)
                     )
    return new{typeof(X), Union{AbsSpec, SpecOpen}, Ring, Hecke.Map}(R)
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
                                      RestrictionType
                                     } <: AbsPreSheaf{
                                                      SpaceType, OpenType,
                                                      OutputType, RestrictionType
                                                     }
  ID::IdDict{AbsSpec, Ideal} # the ideals on the basic patches of the default covering
  OOX::StructureSheafOfRings # the structure sheaf on X
  I::PreSheafOnScheme # the underlying presheaf of ideals for caching

  ### Ideal sheaves on covered schemes
  function IdealSheaf(X::AbsCoveredScheme, ID::IdDict{AbsSpec, Ideal};
      check::Bool=true
    )
    OOX = StructureSheafOfRings(X)

    ### Production of the rings of regular functions; to be cached
    function production_func(F::AbsPreSheaf, U::AbsSpec)
      # If U is an affine chart on which the ideal has already been described, take that.
      haskey(ID, U) && return ID[U]
      # The ideal sheaf has to be provided on at least one dense
      # open subset of every connected component.
      # Otherwise, the ideal sheaf is given by the unit
      # ideals.
      for G in values(glueings(default_covering(space(F))))
        A, B = patches(G)
        Asub, Bsub = glueing_domains(G)
        if A === U && haskey(ID, B) && is_dense(Asub)
          Z = intersect(subscheme(B, ID[B]), Bsub)
          f, _ = glueing_morphisms(G)
          pZ = preimage(f, Z)
          ZU = closure(pZ, U)
          ID[U] = ideal(OO(U), gens(saturated_ideal(modulus(OO(ZU)))))
          return ID[U]
        end
      end
      return ideal(OO(U), one(OO(U)))
    end
    function production_func(F::AbsPreSheaf, U::PrincipalOpenSubset)
      haskey(ID, U) && return ID[U]
      V = ambient_scheme(U)
      IV = F(V)::Ideal
      rho = OOX(V, U)
      IU = ideal(OO(U), rho.(gens(IV)))
      return IU
    end
    function production_func(F::AbsPreSheaf, U::SimplifiedSpec)
      haskey(ID, U) && return ID[U]
      V = original(U)
      IV = F(V)::Ideal
      rho = OOX(V, U)
      IU = ideal(OO(U), rho.(gens(IV)))
      return IU
    end

    ### Production of the restriction maps; to be cached
    function restriction_func(F::AbsPreSheaf, V::AbsSpec, U::AbsSpec)
      return OOX(V, U) # This does not check containment of the arguments
                       # in the ideal. But this is not a parent check and
                       # hence expensive, so we might want to not do that.
    end

    Ipre = PreSheafOnScheme(X, production_func, restriction_func,
                      OpenType=AbsSpec, OutputType=Ideal,
                      RestrictionType=Hecke.Map,
                      is_open_func=_is_open_func_for_schemes_without_specopen(X)
                     )
    I = new{typeof(X), AbsSpec, Ideal, Hecke.Map}(ID, OOX, Ipre)
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
