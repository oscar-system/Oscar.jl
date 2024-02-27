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
  production_func::Function # To produce ℱ(U) for U ⊂ X
  restriction_func::Function # To produce the restriction maps ℱ(U) → ℱ(V) for V ⊂ U ⊂ X open

  function PreSheafOnScheme(X::Scheme, production_func::Any, restriction_func::Any;
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

    ### Production of the rings of regular functions; to be cached
    function production_func(F::AbsPreSheaf, U::AbsAffineScheme)
      return OO(U)
    end
    function production_func(F::AbsPreSheaf, U::AffineSchemeOpenSubscheme)
      return OO(U)
    end

    ### Production of the restriction maps; to be cached
    function restriction_func(F::AbsPreSheaf, V::AbsAffineScheme, U::AbsAffineScheme)
      V === U || error("basic affine patches must be the same")
      return identity_map(OO(V))
    end
    function restriction_func(
        F::AbsPreSheaf, 
        V::AbsAffineScheme, 
        U::Union{<:PrincipalOpenSubset, <:SimplifiedAffineScheme}
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
        f, g = gluing_morphisms(G)
        pbg = pullback(g)
        function rho_func(x::RingElem)
          parent(x) === OV || error("element does not belong to the correct domain")
          y = pbg(domain(pbg)(x, check=false))
          yy = restrict(y, U_flat, check=false)
          return pullback(incU)(yy)
        end
        return hom(OV, OU, rho_func.(gens(OV)), check=false)
      end
      error("arguments are not valid")
    end

    function restriction_func(F::AbsPreSheaf, 
        V::Union{<:PrincipalOpenSubset, <:SimplifiedAffineScheme},
        U::AbsAffineScheme
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
        # kept for the moment because of possible incompatibilities with gluings 
        # along AffineSchemeOpenSubschemes.
        function rho_func(a::RingElem)
          parent(a) === OV || error("element does not belong to the correct ring")
          # We may assume that all denominators admissible in V are
          # already units in OO(U)
          return OU(lifted_numerator(a))*inv(OU(lifted_denominator(a)))
        end
        return hom(OV, OU, rho_func.(gens(OV)), check=false)
      else
        G = default_covering(X)[W, U]
        W1, W2 = gluing_domains(G)
        f, g = gluing_morphisms(G)
        g_res = restrict(g, U, V_direct, check=false)
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
        V::Union{<:PrincipalOpenSubset, <:SimplifiedAffineScheme},
        U::Union{<:PrincipalOpenSubset, <:SimplifiedAffineScheme}
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
                   pullback(inc_U_flat).(pullback(inclusion_morphism(U_flat, V_flat, check=false)).(pullback(inverse(inc_V_flat)).(gens(OV)))), check=false
                  )
      else
        G = default_covering(X)[A, B]
        f, g = gluing_morphisms(G)
        VV_flat = intersect(V_flat, codomain(f))
        VU = preimage(f, VV_flat, check=false)
        fres = restrict(f, VU, VV_flat, check=false)
        inc_V_flat_inv = inverse(inc_V_flat)
        function rho_func(x::RingElem)
          parent(x) === OV || error("input not valid")
          y = pullback(inverse(inc_V_flat))(x)
          y = restrict(y, VV_flat, check=false)
          y = pullback(fres)(y)
          y = restrict(y, U_flat, check=false)
          return pullback(inc_U_flat)(y)
        end
        return hom(OV, OU, rho_func.(gens(OV)), check=false)
      end
      error("arguments are invalid")
    end

    function restriction_func(F::AbsPreSheaf, V::AbsAffineScheme, W::AffineSchemeOpenSubscheme)
      OV = F(V)
      OW = F(W)
      V in default_covering(X) || return false
      ambient_scheme(W) in default_covering(X) || return false
      if V === ambient_scheme(W)
        return MapFromFunc(OV, OW, x->(OW(x)))
      else
        G = default_covering(X)[V, ambient_scheme(W)]
        f, g = gluing_morphisms(G)
        function rho_func(a::RingElem)
          parent(a) === OV || error("element does not belong to the correct ring")
          return restrict(pullback(g)(OO(domain(f))(a)), W, check=false)
        end
        return MapFromFunc(OV, OW, rho_func)
      end
    end

    ### cleaned up until here ###
    # We do not make AffineSchemeOpenSubscheme compatible with the tree structures, yet. 
    # All AffineSchemeOpenSubscheme's are hence required to have an ambient_scheme on the top level. 

    function restriction_func(F::AbsPreSheaf, V::PrincipalOpenSubset, W::AffineSchemeOpenSubscheme)
      error("method not implemented at the moment")
      OV = F(V)
      OW = F(W)
      if ambient_scheme(V) === ambient_scheme(W)
        function rho_func(a::RingElem)
          parent(a) === OV || error("element does not belong to the correct ring")
          return OW(a)
        end
        return MapFromFunc(OV, OW, rho_func)
      else
        G = default_covering(X)(ambient_scheme(V), ambient_scheme(W))
        f, g = gluing_morphisms(G)
        VG = intersect(V, domain(f))
        preV = preimage(g, VG, check=false)
        gres = restriction(g, preV, VG, check=false)
        inc = inclusion_morphism(W, preV, check=false)
        function rho_func2(a::RingElem)
          parent(a) === OV || error("element does not belong to the correct ring")
          return pullback(inc)(pullback(gres)(OO(preV)(a)))
        end
        return MapFromFunc(OV, OW, rho_func2)
      end
    end
    function restriction_func(F::AbsPreSheaf, V::AffineSchemeOpenSubscheme, W::AffineSchemeOpenSubscheme)
      OV = F(V)
      OW = F(W)
      if ambient_scheme(V) === ambient_scheme(W)
        inc = inclusion_morphism(W, V, check=false)
        return MapFromFunc(OV, OW, pullback(inc))
      else
        G = default_covering(X)[ambient_scheme(V), ambient_scheme(W)]
        f, g = gluing_morphisms(G)
        VG = intersect(V, domain(f))
        inc0 = inclusion_morphism(VG, V, check=false)
        preV = preimage(g, VG, check=false)
        gres = restrict(g, preV, VG, check=false)
        inc = inclusion_morphism(W, preV, check=false)
        return MapFromFunc(OV, OW, x->(pullback(inc)(pullback(gres)(pullback(inc0)(x)))))
      end
    end

    R = PreSheafOnScheme(X, production_func, restriction_func,
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
    IdealSheaf <: AbsPreSheaf

A sheaf of ideals ``ℐ`` on an `AbsCoveredScheme` ``X``.

Note that due to technical reasons, the admissible open subsets are restricted
to the following:
 * `U::AbsAffineScheme` among the `basic_patches` of the `default_covering` of `X`;
 * `U::PrincipalOpenSubset` with `ambient_scheme(U)` in the `basic_patches` of the `default_covering` of `X`.

One can call the restriction maps of ``ℐ`` across charts, implicitly using the
identifications given by the gluings in the `default_covering`.
"""
@attributes mutable struct IdealSheaf{SpaceType, OpenType, OutputType,
                                      RestrictionType
                                     } <: AbsPreSheaf{
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
    OOX = StructureSheafOfRings(X)

    function production_func(F::AbsPreSheaf, U::AbsAffineScheme)
      # If U is an affine chart on which the ideal has already been described, take that.
      haskey(ID, U) && return ID[U]
      # Transferring from another chart did not work. That means 
      # I(U) is already prescribed on some refinement of U. We 
      # need to gather that information from all the patches involved
      # and assemble the ideal from there.
      V = [W for W in keys(ID) if has_ancestor(x->(x===U), W)] # gather all patches under U

      # Check for some SimplifiedAffineScheme lurking around
      if any(x->(x isa SimplifiedAffineScheme), V)
        i = findfirst(x->(x isa SimplifiedAffineScheme), V)
        W = V[i]
        _, g = identification_maps(W)
        return ideal(OO(U), pullback(g).(gens(ID[W])))
      end

      length(V) == 0 && return ideal(OO(U), one(OO(U))) # In this case really nothing is defined here.
                                                        # Just return the unit ideal so that the 
                                                        # associated subscheme is empty.
      result = ideal(OO(U), one(OO(U)))
      V = filter!(x->(x isa PrincipalOpenSubset && ambient_scheme(x) === U), V)
      for VV in V
        result = intersect(result, ideal(OO(U), gens(saturated_ideal(production_func(F, VV)))))
      end
      return result
    end
    function production_func(F::AbsPreSheaf, U::PrincipalOpenSubset)
      haskey(ID, U) && return ID[U]
      V = ambient_scheme(U)
      # In case the ambient_scheme is leading out of the admissible domain, 
      # this is a top-chart and we have to reconstruct from below.
      if !is_open_func(F)(V, space(F)) 
        V = [W for W in keys(ID) if has_ancestor(x->(x===U), W)] # gather all patches under U

        # Check for some SimplifiedAffineScheme lurking around
        if any(x->(x isa SimplifiedAffineScheme), V)
          i = findfirst(x->(x isa SimplifiedAffineScheme), V)
          W = V[i]
          _, g = identification_maps(W)
          return ideal(OO(U), pullback(g).(gens(ID[W])))
        end

        length(V) == 0 && return ideal(OO(U), one(OO(U))) # In this case really nothing is defined here.
        # Just return the unit ideal so that the 
        # associated subscheme is empty.
        result = ideal(OO(U), one(OO(U)))
        V = filter!(x->(x isa PrincipalOpenSubset && ambient_scheme(x) === U), V)
        for VV in V
          result = intersect(result, ideal(OO(U), gens(saturated_ideal(production_func(F, VV)))))
        end
        return result
      end

      IV = F(V)::Ideal
      rho = OOX(V, U)
      IU = ideal(OO(U), rho.(gens(IV)))
      return IU
    end
    function production_func(F::AbsPreSheaf, U::SimplifiedAffineScheme)
      haskey(ID, U) && return ID[U]
      V = original(U)
      # In case the original is leading out of the admissible domain, 
      # this is a top-chart and we have to reconstruct from below.
      if !is_open_func(F)(V, space(F)) 
        V = [W for W in keys(ID) if has_ancestor(x->(x===U), W)] # gather all patches under U

        # Check for some SimplifiedAffineScheme lurking around
        if any(x->(x isa SimplifiedAffineScheme), V)
          i = findfirst(x->(x isa SimplifiedAffineScheme), V)
          W = V[i]
          _, g = identification_maps(W)
          return ideal(OO(U), pullback(g).(gens(ID[W])))
        end

        length(V) == 0 && return ideal(OO(U), one(OO(U))) # In this case really nothing is defined here.
        # Just return the unit ideal so that the 
        # associated subscheme is empty.
        result = ideal(OO(U), one(OO(U)))
        V = filter!(x->(x isa PrincipalOpenSubset && ambient_scheme(x) === U), V)
        for VV in V
          result = intersect(result, ideal(OO(U), gens(saturated_ideal(production_func(F, VV)))))
        end
        return result
      end
      IV = F(V)::Ideal
      rho = OOX(V, U)
      IU = ideal(OO(U), rho.(gens(IV)))
      return IU
    end

    ### Production of the restriction maps; to be cached
    function restriction_func(F::AbsPreSheaf, V::AbsAffineScheme, U::AbsAffineScheme)
      return OOX(V, U) # This does not check containment of the arguments
                       # in the ideal. But this is not a parent check and
                       # hence expensive, so we might want to not do that.
    end

    Ipre = PreSheafOnScheme(X, production_func, restriction_func,
                      OpenType=AbsAffineScheme, OutputType=Ideal,
                      RestrictionType=Map,
                      is_open_func=_is_open_func_for_schemes_without_affine_scheme_open_subscheme(X)
                     )
    I = new{typeof(X), AbsAffineScheme, Ideal, Map}(ID, OOX, Ipre)
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
            I(A[i]) == ideal(OOX(A[i]), I(V, A[i]).(gens(I(V)))) || error("ideals do not coincide on overlap")
          end
          for i in 1:number_of_complement_equations(B)
            I(B[i]) == ideal(OOX(B[i]), I(U, B[i]).(gens(I(U)))) || error("ideals do not coincide on overlap")
          end
        end
      end
    end
    return I
  end
end
