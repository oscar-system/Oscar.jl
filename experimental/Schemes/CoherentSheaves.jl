export CoherentSheaf
export tautological_bundle
export twisting_sheaf
export cotangent_sheaf
export tangent_sheaf
export sheaf_of_rings
export dual
export LineBundle
export PushforwardSheaf
export is_locally_free
#export PullbackSheaf

abstract type AbsCoherentSheaf{
                               SpaceType, OpenType,
                               OutputType, RestrictionType
                              } <: AbsPreSheaf{
                                               SpaceType, OpenType,
                                               OutputType, RestrictionType
                                              } 
end

### Interface for coherent sheaves

@Markdown.doc """
    scheme(F::AbsCoherentSheaf)

Return the scheme on which this sheaf is defined.
"""
scheme(F::AbsCoherentSheaf) = space(underlying_presheaf(F))

@Markdown.doc """
    sheaf_of_rings(F::AbsCoherentSheaf) 

Return the sheaf of rings over which ``ℱ`` is defined.
"""
function sheaf_of_rings(F::AbsCoherentSheaf) 
  error("method not implemented for coherent sheaves of type $(typeof(F))")
end

function Base.show(io::IO, M::AbsCoherentSheaf)
  print(io, "sheaf of $(sheaf_of_rings(M))-modules on $(scheme(M))")
end


### The following provides a function for the internal checks that 
# a given set U is open in and admissible for sheaves of modules on X.
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
function _is_open_for_modules(X::AbsCoveredScheme)
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
      gU = preimage(g, U)
      is_subset(gU, V) || return false
    end
    return true
  end
  function is_open_func(U::PrincipalOpenSubset, Y::AbsCoveredScheme)
    return Y === X && ambient_scheme(U) in affine_charts(X)
  end
  function is_open_func(U::AbsSpec, Y::AbsCoveredScheme)
    return Y === X && U in affine_charts(X)
  end
  function is_open_func(Z::AbsCoveredScheme, Y::AbsCoveredScheme)
    return X === Y === Z
  end
  function is_open_func(U::AbsSpec, V::AbsSpec)
    any(x->x===U, affine_charts(X)) || return false
    any(x->x===V, affine_charts(X)) || return false
    G = default_covering(X)[U, V]
    return issubset(U, glueing_domains(G)[1])
  end
  function is_open_func(
      U::AbsSpec,
      V::Union{<:PrincipalOpenSubset, <:SimplifiedSpec}
    )
    issubset(U, V) && return true
    any(x->x===U, affine_charts(X)) || return false
    inc_V_flat = _flatten_open_subscheme(V, default_covering(X))
    A = ambient_scheme(codomain(inc_V_flat))
    Vdirect = codomain(inc_V_flat)
    W = ambient_scheme(Vdirect)
    haskey(glueings(default_covering(X)), (W, U)) || return false # In this case, they are not glued
    G = default_covering(X)[W, U]
    f, g = glueing_morphisms(G)
    pre_V = preimage(g, V)
    return is_subset(U, pre_V)
  end
  function is_open_func(
      U::Union{<:PrincipalOpenSubset, <:SimplifiedSpec}, 
      V::AbsSpec
    )
    any(x->x===V, affine_charts(X)) || return false
    inc_U_flat = _flatten_open_subscheme(U, default_covering(X))
    A = ambient_scheme(codomain(inc_U_flat))
    Udirect = codomain(inc_U_flat)
    W = ambient_scheme(Udirect)
    haskey(glueings(default_covering(X)), (W, V)) || return false # In this case, they are not glued
    G = default_covering(X)[W, V]
    return is_subset(Udirect, glueing_domains(G)[1])
  end
  return is_open_func
end


########################################################################
# Coherent sheaves of modules on covered schemes                       #
########################################################################
@Markdown.doc """
    SheafOfModules <: AbsPreSheaf

A sheaf of modules ``ℳ`` on an `AbsCoveredScheme` ``X``.

Note that due to technical reasons, the admissible open subsets are restricted
to the following:
 * `U::AbsSpec` among the `basic_patches` of the `default_covering` of `X`;
 * `U::PrincipalOpenSubset` with `ambient_scheme(U)` in the `basic_patches` of the `default_covering` of `X`.

One can call the restriction maps of ``ℳ`` across charts implicitly using the
identifications given by the glueings in the `default_covering`.
"""
@attributes mutable struct SheafOfModules{SpaceType, OpenType, 
                                          OutputType,
                                          RestrictionType
                                         } <: AbsCoherentSheaf{
                                                               SpaceType, OpenType,
                                                               OutputType, RestrictionType
                                                              }
  ID::IdDict{AbsSpec, ModuleFP} # the modules on the basic patches of the default covering
  OOX::StructureSheafOfRings # the structure sheaf on X
  M::PreSheafOnScheme # the underlying presheaf of modules for caching

  ### Sheaves of modules on covered schemes
  function SheafOfModules(X::AbsCoveredScheme, 
      MD::IdDict{AbsSpec, ModuleFP}, # A dictionary of modules on the `affine_charts` of `X`
      MG::IdDict{Tuple{AbsSpec, AbsSpec}, MatrixElem}; # A dictionary for pairs `(U, V)` of 
                                                       # `affine_charts` of `X` such that 
                                                       # A = MG[(U, V)] has entries aᵢⱼ with 
                                                       # gᵢ = ∑ⱼ aᵢⱼ ⋅ fⱼ on U ∩ V with gᵢ the 
                                                       # restrictions of the generators of M[U]
                                                       # and fⱼ the restrictions of the generators 
                                                       # of MD[V]. The aᵢⱼ are elements of 𝒪(U ∩ V)
                                                       # represented as a subset of V.
      check::Bool=true
    )
    OOX = OO(X)
    #OOX = StructureSheafOfRings(X)

    ### Production of the modules on open sets; to be cached
    function production_func(
        F::AbsPreSheaf,
        U::AbsSpec
      )
      haskey(MD, U) && return MD[U]
      error("module on $U was not found")
    end

    function production_func(
        F::AbsPreSheaf,
        U::PrincipalOpenSubset
      )
      V = ambient_scheme(U)
      MV = F(V)
      rho = OOX(V, U)
      MU, phi = change_base_ring(rho, MV)
      add_incoming_restriction!(F, V, MU, phi)
      return MU
    end

    function production_func(
        F::AbsPreSheaf,
        U::SimplifiedSpec
      )
      V = original(U)
      MV = F(V)
      rho = OOX(V, U)
      MU, phi = change_base_ring(rho, MV)
      add_incoming_restriction!(F, V, MU, phi)
      return MU
    end

    ### Production of the restriction maps; to be cached
    function restriction_func(F::AbsPreSheaf, V::AbsSpec, U::AbsSpec)
      if any(W->(W === U), affine_charts(X)) && any(W->(W === V), affine_charts(X))
        MV = F(V)
        MU = F(U)
        A = MG[(V, U)] # The transition matrix
        UU, _ = glueing_domains(default_covering(X)[U, V])
        psi = OOX(UU, U) # Needs to exist by the checks of is_open_func, even though 
        # in general UU ⊂ U!
        return hom(MV, MU, [sum([psi(A[i, j]) * MU[j] for j in 1:ngens(MU)], init=zero(MU)) for i in 1:ngens(MV)], OOX(V, U))
      else
        error("invalid input")
      end
    end

    function restriction_func(F::AbsPreSheaf, V::AbsSpec, U::PrincipalOpenSubset)
      # If V was not an affine_chart of X, some other function would have 
      # been triggered. 

      # First the easy case: Inheritance from an ancestor in the tree.
      if ambient_scheme(U) === V
        # If the restriction was more complicated than what follows, then 
        # it would have been cached earlier and this call would not have happened
        # This is the end of the recursion induced in the next elseif below.
        res = hom(F(V), F(U), gens(F(U)), OOX(W, U))
        return res
      elseif some_ancestor(W->(W === V), U)
        W = ambient_scheme(U)
        return compose(F(V, W), F(W, U))
      end

      # Now we know we have a transition across charts
      W = __find_chart(U, default_covering(X))
      A = MG[(V, W)] # The transition matrix
      WW, _ = glueing_domains(default_covering(X)[W, V])
      # From W to U (and hence also from WW to U) the generators of the modules 
      # in F might have changed. Thus, we have to expect a non-trivial transition 
      # from the top-level down to U. The transition matrix A is only given with 
      # respect to the generators of F(W), so we have to map them manually down.
      # The call to F(W, U) will be handled by the above if-clauses.
      return hom(F(V), F(U), 
                 [sum([OOX(WW, U)(A[i, j])*F(W, U)(F(W)[j]) for j in 1:ngens(F(W))], init=zero(F(U))) 
                  for i in 1:ngens(F(V))], 
                 OOX(V, U)
                )
    end

    function restriction_func(F::AbsPreSheaf, V::AbsSpec, U::SimplifiedSpec)
      # If V was not an affine_chart of X, some other function would have 
      # been triggered. 

      # First the easy case: Inheritance from an ancestor in the tree.
      if original(U) === V
        # If the restriction was more complicated than what follows, then 
        # it would have been cached earlier and this call would not have happened
        # This is the end of the recursion induced in the next elseif below.
        W = original(U)
        res = hom(F(W), F(U), gens(MU), OOX(W, U))
        return res
      elseif some_ancestor(W->(W === V), U)
        W = original(U)
        return compose(F(V, W), F(W, U))
      end

      # Now we know we have a transition across charts
      W = __find_chart(U, default_covering(X))
      A = MG[(V, W)] # The transition matrix
      WW, _ = glueing_domains(default_covering(X)[W, V])
      # From W to U (and hence also from WW to U) the generators of the modules 
      # in F might have changed. Thus, we have to expect a non-trivial transition 
      # from the top-level down to U. The transition matrix A is only given with 
      # respect to the generators of F(W), so we have to map them manually down.
      # The call to F(W, U) will be handled by the above if-clauses.
      return hom(F(V), F(U), 
                 [sum([OOX(WW, U)(A[i, j])*F(W, U)(F(W)[j]) for j in 1:ngens(F(W))], init=zero(F(U)))
                  for i in 1:ngens(F(V))], 
                 OOX(V, U)
                )
    end
    function restriction_func(F::AbsPreSheaf, 
        V::Union{<:PrincipalOpenSubset, <:SimplifiedSpec}, 
        U::AbsSpec
      )
      # We know that V can not be an ancestor of U, but U must be an affine chart.
      # Probably even an ancestor of V itself. 
      W = __find_chart(V, default_covering(X))
      if W === U
        # U and V must actually be isomorphic, but the modules of F might be 
        # represented in different ways. We have to construct the inverse of 
        # the restriction map from U to V.
        gens_U = F(U, V).(gens(F(U)))
        M, inc = sub(F(V), gens_U)
        img_gens = elem_type(F(U))[]
        for v in gens(F(V))
          w = preimage(inc, v)
          c = coordinates(w)
          push!(img_gens,
                sum(OOX(V, U)(c[i])*gens(F(U), i) for i in 1:ngens(F(U)))
               )
        end
        return hom(F(V), F(U), img_gens, OOX(V, U))
      else
        # U must be properly contained in the glueing domains of the 
        # glueing of the affine chart of V with U.
        error("case not implemented")
      end
      # Problem: We can assume that we know how to pass from generators 
      # of W = __find_chart(V, default_covering(X)) to those on V, but we do not 
      # know the inverse to this. But the transition matrix to U is given 
      # with respect to the generators on W.
      error("case not implemented")
    end
    function restriction_func(F::AbsPreSheaf, V::PrincipalOpenSubset, U::PrincipalOpenSubset)
      V === U && return identity_map(F(U))

      if V === ambient_scheme(U)
        return hom(F(V), F(U), gens(F(U)), OOX(V, U)) # If this had been more complicated, it would have been cached.
      elseif some_ancestor(W->W===V, U)
        W = ambient_scheme(U)
        return compose(F(V, W), F(W, U))
      end

      # Below follow the more complicated cases. 
      success, _ = _have_common_ancestor(U, V)
      if success
        W = __find_chart(U, default_covering(X))
        gens_U = F(W, U).(gens(F(W))) # This will be caught by the preceeding clauses
        gens_V = F(W, V).(gens(F(W)))
        sub_V, inc = sub(F(V), gens_V)
        img_gens = elem_type(F(U))[]
        for v in gens(F(V))
          w = preimage(inc, v) # We know that inc is actually an isomorphism
          c = coordinates(w)
          w = sum(OOX(V, U)(c[i])*gens_U[i] 
                  for i in 1:length(gens_U)
                 )
          push!(img_gens, w)
        end
        return hom(F(V), F(U), img_gens, OOX(V, U))
      end

      # Now we know we have a transition between different charts.
      inc_U = _flatten_open_subscheme(U, default_covering(X))
      inc_V = _flatten_open_subscheme(V, default_covering(X))
      U_flat = codomain(inc_U)
      V_flat = codomain(inc_V)
      WU = ambient_scheme(U_flat)
      WV = ambient_scheme(V_flat)
      WU = __find_chart(U, default_covering(X))
      WV = __find_chart(V, default_covering(X))
      # The problem is: The generators of F(WU) may be different from 
      # those of F(U) and similarly for V. But the transition matrices 
      # are only described for those on WU and WV. Thus we need to 
      # implicitly do a base change. This is done by forwarding the generators 
      # of F(WU) to F(U) and expressing it in terms of the generators there. 
      gens_U = F(WU, U).(gens(F(WU))) # This will be caught by the preceeding clauses
      gens_V = F(WV, V).(gens(F(WV)))
      sub_V, inc = sub(F(V), gens_V)
      img_gens = elem_type(F(U))[]
      A = MG[(WV, WU)] # The transition matrix
      WW, _ = glueing_domains(default_covering(X)[WU, WV])
      for v in gens(F(V))
        w = preimage(inc, v) # We know that inc is actually an isomorphism
        c = coordinates(w)
        w = sum(sum(OOX(V, U)(c[i])*OOX(WW, U)(A[i, j])*gens_U[j] 
                    for i in 1:length(gens_V)) 
                for j in 1:length(gens_U)
               )
        push!(img_gens, w)
      end
      return hom(F(V), F(U), img_gens, OOX(V, U))
    end
    function restriction_func(F::AbsPreSheaf, V::PrincipalOpenSubset, U::SimplifiedSpec)
      if V === original(U)
        return hom(F(V), F(U), gens(F(U)), OOX(V, U)) # If this had been more complicated, it would have been cached.
      elseif some_ancestor(W->W===V, U)
        W = original(U)
        return compose(F(V, W), F(W, U))
      end

      # Below follow the more complicated cases. 
      success, _ = _have_common_ancestor(U, V)
      if success
        W = __find_chart(U, default_covering(X))
        gens_U = F(W, U).(gens(F(W))) # This will be caught by the preceeding clauses
        gens_V = F(W, V).(gens(F(W)))
        sub_V, inc = sub(F(V), gens_V)
        img_gens = elem_type(F(U))[]
        for v in gens(F(V))
          w = preimage(inc, v) # We know that inc is actually an isomorphism
          c = coordinates(w)
          w = sum(OOX(V, U)(c[i])*gens_U[i] 
                  for i in 1:length(gens_U)
                 )
          push!(img_gens, w)
        end
        return hom(F(V), F(U), img_gens, OOX(V, U))
      end

      # Now we know we have a transition between different charts.
      inc_U = _flatten_open_subscheme(U, default_covering(X))
      inc_V = _flatten_open_subscheme(V, default_covering(X))
      U_flat = codomain(inc_U)
      V_flat = codomain(inc_V)
      WU = ambient_scheme(U_flat)
      WV = ambient_scheme(V_flat)
      WU = __find_chart(U, default_covering(X))
      WV = __find_chart(V, default_covering(X))
      # The problem is: The generators of F(WU) may be different from 
      # those of F(U) and similarly for V. But the transition matrices 
      # are only described for those on WU and WV. Thus we need to 
      # implicitly do a base change. This is done by forwarding the generators 
      # of F(WU) to F(U) and expressing it in terms of the generators there. 
      gens_U = F(WU, U).(gens(F(WU))) # This will be caught by the preceeding clauses
      gens_V = F(WV, V).(gens(F(WV)))
      sub_V, inc = sub(F(V), gens_V)
      img_gens = elem_type(F(U))[]
      A = MG[(WV, WU)] # The transition matrix
      WW, _ = glueing_domains(default_covering(X)[WU, WV])
      for v in gens(F(V))
        w = preimage(inc, v) # We know that inc is actually an isomorphism
        c = coordinates(w)
        w = sum(sum(OOX(V, U)(c[i])*OOX(WW, U)(A[i, j])*gens_U[j] 
                    for i in 1:length(gens_V)) 
                for j in 1:length(gens_U)
               )
        push!(img_gens, w)
      end
      return hom(F(V), F(U), img_gens, OOX(V, U))
    end
    function restriction_func(F::AbsPreSheaf, V::SimplifiedSpec, U::PrincipalOpenSubset)
      if V === ambient_scheme(U)
        return hom(F(V), F(U), gens(F(U)), OOX(V, U)) # If this had been more complicated, it would have been cached.
      elseif some_ancestor(W->W===V, U)
        W = ambient_scheme(U)
        return compose(F(V, W), F(W, U))
      end

      # Below follow the more complicated cases. 
      success, _ = _have_common_ancestor(U, V)
      if success
        W = __find_chart(U, default_covering(X))
        gens_U = F(W, U).(gens(F(W))) # This will be caught by the preceeding clauses
        gens_V = F(W, V).(gens(F(W)))
        sub_V, inc = sub(F(V), gens_V)
        img_gens = elem_type(F(U))[]
        for v in gens(F(V))
          w = preimage(inc, v) # We know that inc is actually an isomorphism
          c = coordinates(w)
          w = sum(OOX(V, U)(c[i])*gens_U[i] 
                  for i in 1:length(gens_U)
                 )
          push!(img_gens, w)
        end
        return hom(F(V), F(U), img_gens, OOX(V, U))
      end

      # Now we know we have a transition between different charts.
      inc_U = _flatten_open_subscheme(U, default_covering(X))
      inc_V = _flatten_open_subscheme(V, default_covering(X))
      U_flat = codomain(inc_U)
      V_flat = codomain(inc_V)
      WU = ambient_scheme(U_flat)
      WV = ambient_scheme(V_flat)
      WU = __find_chart(U, default_covering(X))
      WV = __find_chart(V, default_covering(X))
      # The problem is: The generators of F(WU) may be different from 
      # those of F(U) and similarly for V. But the transition matrices 
      # are only described for those on WU and WV. Thus we need to 
      # implicitly do a base change. This is done by forwarding the generators 
      # of F(WU) to F(U) and expressing it in terms of the generators there. 
      gens_U = F(WU, U).(gens(F(WU))) # This will be caught by the preceeding clauses
      gens_V = F(WV, V).(gens(F(WV)))
      sub_V, inc = sub(F(V), gens_V)
      img_gens = elem_type(F(U))[]
      A = MG[(WV, WU)] # The transition matrix
      WW, _ = glueing_domains(default_covering(X)[WU, WV])
      for v in gens(F(V))
        w = preimage(inc, v) # We know that inc is actually an isomorphism
        c = coordinates(w)
        w = sum(sum(OOX(V, U)(c[i])*OOX(WW, U)(A[i, j])*gens_U[j] 
                    for i in 1:length(gens_V)) 
                for j in 1:length(gens_U)
               )
        push!(img_gens, w)
      end
      return hom(F(V), F(U), img_gens, OOX(V, U))
    end
    function restriction_func(F::AbsPreSheaf, V::SimplifiedSpec, U::SimplifiedSpec)
      V === U && return identity_map(F(U))

      if V === original(U)
        return hom(F(V), F(U), gens(F(U)), OOX(V, U)) # If this had been more complicated, it would have been cached.
      elseif some_ancestor(W->W===V, U)
        W = original(U)
        return compose(F(V, W), F(W, U))
      end

      # Below follow the more complicated cases. 
      success, _ = _have_common_ancestor(U, V)
      if success
        W = __find_chart(U, default_covering(X))
        gens_U = F(W, U).(gens(F(W))) # This will be caught by the preceeding clauses
        gens_V = F(W, V).(gens(F(W)))
        sub_V, inc = sub(F(V), gens_V)
        img_gens = elem_type(F(U))[]
        for v in gens(F(V))
          w = preimage(inc, v) # We know that inc is actually an isomorphism
          c = coordinates(w)
          w = sum(OOX(V, U)(c[i])*gens_U[i] 
                  for i in 1:length(gens_U)
                 )
          push!(img_gens, w)
        end
        return hom(F(V), F(U), img_gens, OOX(V, U))
      end

      # Now we know we have a transition between different charts.
      inc_U = _flatten_open_subscheme(U, default_covering(X))
      inc_V = _flatten_open_subscheme(V, default_covering(X))
      U_flat = codomain(inc_U)
      V_flat = codomain(inc_V)
      WU = ambient_scheme(U_flat)
      WV = ambient_scheme(V_flat)
      WU = __find_chart(U, default_covering(X))
      WV = __find_chart(V, default_covering(X))
      # The problem is: The generators of F(WU) may be different from 
      # those of F(U) and similarly for V. But the transition matrices 
      # are only described for those on WU and WV. Thus we need to 
      # implicitly do a base change. This is done by forwarding the generators 
      # of F(WU) to F(U) and expressing it in terms of the generators there. 
      gens_U = F(WU, U).(gens(F(WU))) # This will be caught by the preceeding clauses
      gens_V = F(WV, V).(gens(F(WV)))
      sub_V, inc = sub(F(V), gens_V)
      img_gens = elem_type(F(U))[]
      A = MG[(WV, WU)] # The transition matrix
      WW, _ = glueing_domains(default_covering(X)[WU, WV])
      for v in gens(F(V))
        w = preimage(inc, v) # We know that inc is actually an isomorphism
        c = coordinates(w)
        w = sum(sum(OOX(V, U)(c[i])*OOX(WW, U)(A[i, j])*gens_U[j] 
                    for i in 1:length(gens_V)) 
                for j in 1:length(gens_U)
               )
        push!(img_gens, w)
      end
      return hom(F(V), F(U), img_gens, OOX(V, U))
    end

    Mpre = PreSheafOnScheme(X, production_func, restriction_func,
                      OpenType=AbsSpec, OutputType=ModuleFP,
                      RestrictionType=Hecke.Map,
                      is_open_func=_is_open_func_for_schemes_without_specopen(X)
                      #is_open_func=_is_open_for_modules(X)
                     )
    M = new{typeof(X), AbsSpec, ModuleFP, Hecke.Map}(MD, OOX, Mpre)
    if check
      # Check that all sheaves of modules are compatible on the overlaps.
      # TODO: eventually replace by a check that on every basic
      # affine patch, the ideal sheaf can be inferred from what is
      # given on one dense open subset.
    end
    return M
  end
end

### forwarding and implementing the required getters
underlying_presheaf(M::SheafOfModules) = M.M
sheaf_of_rings(M::SheafOfModules) = M.OOX


@Markdown.doc """
    twisting_sheaf(IP::ProjectiveScheme{<:Field}, d::Int)

For a `ProjectiveScheme` ``ℙ`` return the ``d``-th twisting sheaf 
``𝒪(d)`` as a `CoherentSheaf` on ``ℙ``.
"""
function twisting_sheaf(IP::ProjectiveScheme{<:Field}, d::Int)
  # First, look up whether this sheaf has already been computed:
  if !has_attribute(IP, :twisting_sheaves)
    set_attribute!(IP, :twisting_sheaves, Dict{Int, SheafOfModules}())
  end
  twisting_sheaves = get_attribute(IP, :twisting_sheaves)
  haskey(twisting_sheaves, d) && return twisting_sheaves[d]

  X = covered_scheme(IP)
  MD = IdDict{AbsSpec, ModuleFP}()
  S = ambient_coordinate_ring(IP)
  n = ngens(S)-1
  for i in 1:n+1
    U = affine_charts(X)[i]
    MD[U] = FreeMod(OO(U), ["$(symbols(S)[i])^$d"])
  end

  MG = IdDict{Tuple{AbsSpec, AbsSpec}, MatrixElem}()
  C = default_covering(X)
  for G in values(glueings(C))
    (U, V) = patches(G)
    (UU, VV) = glueing_domains(G)
    h_U = complement_equation(UU)
    h_V = complement_equation(VV)
    MG[(U, V)] = (d>= 0 ? OO(VV)(h_V^d) : (inv(OO(VV)(h_V))^(-d)))*one(MatrixSpace(OO(VV), 1, 1))
    MG[(V, U)] = (d>= 0 ? OO(UU)(h_U^d) : (inv(OO(UU)(h_U))^(-d)))*one(MatrixSpace(OO(UU), 1, 1))
  end

  M = SheafOfModules(X, MD, MG)
  # Cache the result for the next usage
  twisting_sheaves[d] = M
  return M
end

@Markdown.doc """
    tautological_bundle(IP::ProjectiveScheme{<:Field})

For a `ProjectiveScheme` ``ℙ`` return the sheaf ``𝒪(-1)`` as a `CoherentSheaf` on ``ℙ``.
"""
function tautological_bundle(IP::ProjectiveScheme{<:Field})
    return twisting_sheaf(IP, -1)
end

@Markdown.doc """
    cotangent_sheaf(X::AbsCoveredScheme)

For an `AbsCoveredScheme` ``X``, return the sheaf ``Ω¹(X)`` of Kaehler-differentials 
on ``X`` as a `CoherentSheaf`.

"""
@attr SheafOfModules function cotangent_sheaf(X::AbsCoveredScheme)
  MD = IdDict{AbsSpec, ModuleFP}()
  for U in affine_charts(X)
    MD[U] = cotangent_module(U)
  end
  MG = IdDict{Tuple{AbsSpec, AbsSpec}, MatrixElem}()
  C = default_covering(X)
  for G in values(glueings(C))
    (U, V) = patches(G)
    (UU, VV) = glueing_domains(G)
    (f, g) = glueing_morphisms(G)
    MG[(U, V)] = transpose(jacobi_matrix(pullback(g).(gens(OO(UU)))))
    MG[(V, U)] = transpose(jacobi_matrix(pullback(f).(gens(OO(VV)))))
  end


  M = SheafOfModules(X, MD, MG)
  return M
end

@Markdown.doc """
    cotangent_module(X::AbsSpec)

Return the ``𝒪(X)``-module ``Ω¹(X)`` of Kaehler-differentials on ``X``.
"""
function cotangent_module(X::AbsSpec)
  error("method not implemented for this type of ring")
end

@attr ModuleFP function cotangent_module(X::AbsSpec{<:Field, <:MPolyRing})
  R = OO(X)
  F = FreeMod(R, ["d$(x)" for x in symbols(R)])
  return F
end

@attr ModuleFP function cotangent_module(X::AbsSpec{<:Field, <:MPolyLocalizedRing})
  R = OO(X)
  P = base_ring(R)
  F = FreeMod(R, ["d$(x)" for x in symbols(P)])
  return F
end

@attr ModuleFP function cotangent_module(X::AbsSpec{<:Field, <:MPolyQuo})
  R = OO(X)
  P = base_ring(R)
  F = FreeMod(R, ["d$(x)" for x in symbols(P)])
  rels, _ = sub(F, transpose(change_base_ring(R, jacobi_matrix(gens(modulus(R))))))
  M, _ = quo(F, rels)
  return M
end

@attr ModuleFP function cotangent_module(X::AbsSpec{<:Field, <:MPolyQuoLocalizedRing})
  R = OO(X)
  P = base_ring(R)
  F = FreeMod(R, ["d$(x)" for x in symbols(P)])
  rels, _ = sub(F, transpose(change_base_ring(R, jacobi_matrix(gens(modulus(R))))))
  M, _ = quo(F, rels)
  return M
end

########################################################################
# Hom-Sheaves                                                          #
########################################################################
#= 
  Hom sheaves ℋom(ℱ, 𝒢) are special. 

  First of all, they can be made completely lazy, as their modules 
  on U ⊂ X can be created from ℱ(U) and 𝒢(U) on the fly and the same
  holds for their transition- and restriction functions.

  Second, Hom sheaves come with a domain, a codomain, and an 
  interpretation mapping of their sections as homomorphisms.

  We realize hom sheaves in this way, taking only ℱ and 𝒢 as input 
  in the constructor.
=#

@attributes mutable struct HomSheaf{SpaceType, OpenType, OutputType,
                                    RestrictionType
                                   } <: AbsCoherentSheaf{
                                                         SpaceType, OpenType,
                                                         OutputType, RestrictionType
                                                        }
  domain::AbsCoherentSheaf{SpaceType, OpenType, OutputType, RestrictionType}
  codomain::AbsCoherentSheaf{SpaceType, OpenType, OutputType, RestrictionType}
  OOX::StructureSheafOfRings
  M::PreSheafOnScheme

  function HomSheaf(F::AbsCoherentSheaf, G::AbsCoherentSheaf)
    X = scheme(F)
    X === scheme(G) || error("sheaves must be defined over the same scheme")
    OOX = sheaf_of_rings(F)
    OOX === sheaf_of_rings(G) || error("sheaves must be defined over the same sheaves of rings")

    ### Production of the modules on open sets; to be cached
    function production_func(FF::AbsPreSheaf, U::AbsSpec)
      return hom(F(U), G(U))[1]
    end

    function restriction_func(FF::AbsPreSheaf, V::AbsSpec, U::AbsSpec)
      MV = FF(V)
      MU = FF(U)
      dom_res = F(V, U)
      cod_res = G(V, U)
      f = gens(F(V))
      rf = dom_res.(f)
      # The following two lines will work, because a set of generators for ℱ(V) 
      # always restricts to a set of generators for ℱ(U). Due to changes of 
      # charts, this might be a non-trivial change of bases, however.
      dom_sub, inc = sub(F(U), rf)
      B = [coordinates(e, dom_sub) for e in ambient_representatives_generators(F(U))]
      images = elem_type(MU)[]
      for phi in gens(MV)
        phi_map = element_to_homomorphism(phi)
        images_f = [sum([B[i][j]*cod_res(phi_map(f[j])) for j in 1:length(f)], init=zero(G(U))) for i in 1:length(B)]
        psi = hom(F(U), G(U), images_f)
        push!(images, homomorphism_to_element(MU, psi))
      end

      return hom(MV, MU, images, OOX(V, U)) # TODO: Set check=false?
    end
      
    Mpre = PreSheafOnScheme(X, production_func, restriction_func,
                      OpenType=AbsSpec, OutputType=ModuleFP,
                      RestrictionType=Hecke.Map,
                      is_open_func=_is_open_func_for_schemes_without_specopen(X)
                     )
    M = new{typeof(X), AbsSpec, ModuleFP, Hecke.Map}(F, G, OOX, Mpre)

    return M
  end
end

### forwarding and implementation of the essential getters
underlying_presheaf(M::HomSheaf) = M.M
sheaf_of_rings(M::HomSheaf) = M.OOX
domain(M::HomSheaf) = M.domain
codomain(M::HomSheaf) = M.codomain

function Base.show(io::IO, M::HomSheaf)
  print(io, "sheaf of homomorphisms from $(domain(M)) to $(codomain(M))")
end

@Markdown.doc """
    free_module(R::StructureSheafOfRings, n::Int)

Return the sheaf of free ``𝒪``-modules ``𝒪ⁿ`` for a structure 
sheaf of rings ``𝒪 = R``.
"""
function free_module(R::StructureSheafOfRings, n::Int)
  return free_module(R, ["e_$i" for i in 1:n])
end

function free_module(R::StructureSheafOfRings, gen_names::Vector{String})
  return free_module(R, Symbol.(gen_names))
end

function free_module(R::StructureSheafOfRings, gen_names::Vector{Symbol})
  X = space(R)
  n = length(gen_names)
  MD = IdDict{AbsSpec, ModuleFP}()
  for U in affine_charts(X)
    MD[U] = FreeMod(OO(U), gen_names)
  end

  MG = IdDict{Tuple{AbsSpec, AbsSpec}, MatrixElem}()
  C = default_covering(X)
  for G in values(glueings(C))
    (U, V) = patches(G)
    (UU, VV) = glueing_domains(G)
    MG[(U, V)] = one(MatrixSpace(OO(VV), n, n))
    MG[(V, U)] = one(MatrixSpace(OO(UU), n, n))
  end

  M = SheafOfModules(X, MD, MG)
  return M
end

@Markdown.doc """
    dual(M::SheafOfModules)

For a `SheafOfModules` ``ℳ`` on an `AbsCoveredScheme` ``X``, return 
the ``𝒪_X``-dual ``ℋ om_{𝒪_X}(ℳ , 𝒪_X)`` of ``ℳ``.
"""
@attr function dual(M::SheafOfModules)
  OOX = sheaf_of_rings(M)
  F = free_module(OOX, ["1"])
  return HomSheaf(M, F)
end

@Markdown.doc """
    tangent_sheaf(X::AbsCoveredScheme)

Return the tangent sheaf ``T_X`` of `X`, constructed as ``ℋ om_{𝒪_X}(Ω¹_X, 𝒪_X)``.
"""
@attr HomSheaf function tangent_sheaf(X::AbsCoveredScheme)
  return dual(cotangent_sheaf(X))
end

########################################################################
# Pushforwards of sheaves for closed embeddings                        #
########################################################################

#=
# Let f : X ↪ Y be a closed embedding with ideal sheaf ℐ on Y and ℳ 
# a sheaf of modules on X. For an open set U ⊂ Y we have 
# f_* ℳ (U) to be the 𝒪_X(f⁻¹(U))-module ℳ (f⁻¹(U)), but seen as 
# an 𝒪_Y(U)-module via the natural restriction of functions. 
#
# Mathematically, this is almost an implicit operation. But since we 
# do not have natural bi-module structures, we need to set up a new 
# sheaf of modules on Y, together with canonical identifications 
# with the modules on X. 
#
# It is clear that this can and should be made lazy.
#                                                                     =#

@attributes mutable struct PushforwardSheaf{SpaceType, OpenType, OutputType,
                                            RestrictionType
                                           } <: AbsCoherentSheaf{
                                                                 SpaceType, OpenType,
                                                                 OutputType, RestrictionType
                                                                }
  inc::CoveredClosedEmbedding
  OOX::StructureSheafOfRings
  OOY::StructureSheafOfRings
  M::AbsCoherentSheaf
  ident::IdDict{AbsSpec, Union{Hecke.Map, Nothing}} # a dictionary caching the identifications
  F::PreSheafOnScheme

  function PushforwardSheaf(inc::CoveredClosedEmbedding, M::AbsCoherentSheaf)
    X = domain(inc)
    X === scheme(M) || error("sheaf must be defined over the domain of the embedding")
    OOX = sheaf_of_rings(M)
    Y = codomain(inc)
    OOY = OO(Y)

    ### Production of the modules on open sets; to be cached
    function production_func(FF::AbsPreSheaf, U::AbsSpec)
      # In case X was empty, return the zero module and store nothing in the identifications.
      if isempty(X) 
        ident[U] = nothing
        return FreeMod(OOY(U), 0)
      end

      # Check whether U ⊂ Y has a nontrivial preimage in X
      f = maps_with_given_codomain(inc, U) # there should be at most one!
      if iszero(length(f))
        ident[U] = nothing
        return FreeMod(OOY(U), 0)
      end
      ff = first(f)
      UX = domain(ff)
      MU, ident_map = _pushforward(pullback(ff), image_ideal(ff), M(UX))
      ident[U] = ident_map
      return MU
    end
    function production_func(FF::AbsPreSheaf, U::PrincipalOpenSubset)
      # In case X was empty, return the zero module and store nothing in the identifications.
      if isempty(X) 
        ident[U] = nothing
        return FreeMod(OOY(U), 0)
      end

      # Check whether U ⊂ Y has a nontrivial preimage in X
      f = maps_with_given_codomain(inc, U) # there should be at most one!
      if !iszero(length(f))
        # in this case, we can produce directly from the source
        ff = first(f)
        UX = domain(ff)
        MU, ident_map = _pushforward(pullback(ff), image_ideal(ff), M(UX))
        ident[U] = ident_map
        return MU
      end

      # We need to restrict from the parent
      W = ambient_scheme(U)
      MW = FF(W)
      MU, res = change_base_ring(OOY(W, U), MW)
      add_incoming_restriction!(FF, W, MU, res)
      return MU
    end


    function restriction_func(FF::AbsPreSheaf, V::AbsSpec, U::AbsSpec)
      MYV = FF(V)
      MYU = FF(U)
      incV_list = maps_with_given_codomain(inc, V)
      incU_list = maps_with_given_codomain(inc, U)
      # return the zero homomorphism in case one of the two sets has 
      # empty preimage.
      if iszero(length(incV_list)) || iszero(length(incU_list)) 
        return hom(MYV, MYU, elem_type(MYU)[zero(MYU) for i in 1:ngens(MYV)], OOY(V, U))
      end
      incV = first(incV_list)
      incU = first(incU_list)
      res_orig = M(domain(incV), domain(incU))
      img_gens = res_orig.(gens(M(domain(incV))))
      return hom(MYV, MYU, (x->preimage(ident[U], x)).(img_gens), OOY(V, U))
    end
    
    ident = IdDict{AbsSpec, Union{Hecke.Map, Nothing}}()

    Blubber = PreSheafOnScheme(Y, production_func, restriction_func,
                      OpenType=AbsSpec, OutputType=ModuleFP,
                      RestrictionType=Hecke.Map,
                      is_open_func=_is_open_func_for_schemes_without_specopen(Y)
                      #is_open_func=_is_open_for_modules(Y)
                     )
    MY = new{typeof(Y), AbsSpec, ModuleFP, Hecke.Map}(inc, OOX, OOY, M, ident, Blubber)
    return MY
  end
end

### forwarding and implementing the required getters
underlying_presheaf(M::PushforwardSheaf) = M.F
sheaf_of_rings(M::PushforwardSheaf) = M.OOY
original_sheaf(M::PushforwardSheaf) = M.M
map(M::PushforwardSheaf) = M.inc

function Base.show(io::IO, M::PushforwardSheaf)
  print(io, "pushforward of $(original_sheaf(M)) along $(map(M))")
end

########################################################################
# pushforward of modules                                               #
########################################################################
#
# It is assumed that f : R → S is a map of rings such that S ≅ R/I and 
# M is an S-module. We transform M into an R-module by adding the 
# necessary relations. The return value is that new module M', together 
# with its identification map M' → M. Note that we can not give the 
# inverse of this map, since there is no well-defined underlying ring 
# homomorphism.
function _pushforward(f::Hecke.Map{<:Ring, <:Ring}, I::Ideal, M::FreeMod)
  R = domain(f)
  S = codomain(f)
  base_ring(I) === R || error("ideal is not defined over the correct ring")
  base_ring(M) === S || error("module is not defined over the correct ring")
  FR = FreeMod(R, M.S) # M.S are the symbols of M
  #FRtoM = hom(FR, M, gens(M), f)
  MR, res = quo(FR, (I*FR)[1])
  ident = hom(MR, M, gens(M), f)
  return MR, ident
end

function _pushforward(f::Hecke.Map{<:Ring, <:Ring}, I::Ideal, M::SubQuo)
  R = domain(f)
  S = codomain(f)
  base_ring(I) === R || error("ideal is not defined over the correct ring")
  base_ring(M) === S || error("module is not defined over the correct ring")
  FS = ambient_free_module(M)
  Mgens = ambient_representatives_generators(M)
  Mrels = relations(M)
  FR, identF = _pushforward(f, I, FS)
  MRgens = [preimage(identF, v) for v in ambient_representatives_generators(M)]
  MRrels = elem_type(FR)[preimage(identF, w) for w in relations(M)]
  MRrels_ext = vcat(MRrels, elem_type(FR)[g*e for g in gens(I) for e in gens(FR)])
  MR = quo(sub(FR, MRgens)[1], sub(FR, MRrels_ext)[1])[1]
  ident = hom(MR, M, gens(M), f)
  return MR, ident
end

@attr Bool function is_locally_free(M::AbsCoherentSheaf)
  return all(U->is_projective(M(U))[1], affine_charts(scheme(M)))
end

#@attr Covering function trivializing_covering(M::AbsCoherentSheaf)
@attr function trivializing_covering(M::AbsCoherentSheaf)
  X = scheme(M)
  OOX = OO(X)
  patch_list = Vector{AbsSpec}()
  for U in affine_charts(X)
    patch_list = vcat(patch_list, _trivializing_covering(M, U))
  end
  C = Covering(patch_list)
  fill_with_lazy_glueings!(C, X)
  return C
end

@attr function trivializing_covering(M::HomSheaf)
  X = scheme(M)
  OOX = OO(X)
  # The problem is that every module of a HomSheaf must know that it is 
  # a hom-module. Hence, the way to go is to pass through a common 
  # refinement of domain and codomain and recreate all the hom modules 
  # as free modules on this covering. 
  #
  # But for this, it is not yet clear where to locate the patches of these 
  # refinements in the tree and how to deal with the restriction maps in a 
  # clean way. Say M = Hom(F, G) where F is trivialized on {Uᵢ} and G on 
  # {Vⱼ}. Then W = Uᵢ∩ Vⱼ would have to be a PrincipalOpenSubset of both 
  # Uᵢ and Vⱼ for the restrictions of F and G to induce the proper job 
  # on restrictions to M(W) automatically. 
  # Hence, we need to manually prescribe how to trivialize and restrict 
  # M on the Ws.
  dom_triv = trivializing_covering(domain(M))
  cod_triv = trivializing_covering(codomain(M))
  patch_list = AbsSpec[]
  for U in patches(dom_triv)
    for V in patches(cod_triv)
      success, W = _have_common_ancestor(U, V)
      if success
        incU = _flatten_open_subscheme(U, W)
        incV = _flatten_open_subscheme(V, W)
        UV = intersect(codomain(incU), codomain(incV))::PrincipalOpenSubset
        push!(patch_list, UV)

        dom_UV, dom_res = change_base_ring(OOX(U, UV), domain(M)(U))
        add_incoming_restriction!(domain(M), U, dom_UV, dom_res)
        add_incoming_restriction!(domain(M), W, dom_UV, 
                                  compose(domain(M)(W, U), dom_res))
        object_cache(domain(M))[UV] = dom_UV

        cod_UV, cod_res = change_base_ring(OOX(V, UV), codomain(M)(V))
        add_incoming_restriction!(codomain(M), V, cod_UV, cod_res)
        add_incoming_restriction!(codomain(M), W, cod_UV, 
                                  compose(codomain(M)(W, V), cod_res))
        object_cache(codomain(M))[UV] = cod_UV

        MUV = M(UV) # This will be a free module; we need to prescribe the restrictions!
        MW = M(W)
        img_gens = elem_type(MUV)[]
        # every generator g of MW is a homomorphism. It takes an element 
        # v ∈ domain(M)(W) to w = ϕ_{g}(v) ∈ codomain(M)(W). 
        # Where does g map to when restricting to MUV? 
        #
        for g in gens(MW)
          phi = element_to_homomorphism(g)
          img_gens_phi = cod_res.(codomain(M)(W, V).(phi.(gens(domain(M)(W)))))
          sub_dom, inc_dom = sub(domain(M)(UV), 
                                 domain(M)(W, UV).(gens(domain(M)(W))))
                                 #dom_res.(domain(M)(W, U).(gens(domain(M)(W)))))
          img_gens_psi = elem_type(codomain(M)(UV))[]
          for v in gens(domain(M)(UV))
            w = preimage(inc_dom, v)
            c = coordinates(w) # These are the coordinates in the original set 
                               # of generators in the domain
            # We use this to compute the image of v
            phi_v = sum([c[i]*img_gens_phi[i] for i in 1:length(img_gens_phi)], init=zero(codomain(M)(UV)))
            # and push it to the list.
            push!(img_gens_psi, phi_v)
          end
          # From that list, we can assemble what the restriction of phi 
          # looks like as a homomorphism 
          psi = hom(domain(M)(UV), codomain(M)(UV), img_gens_psi)
          # and convert it to a module element.
          img_g = homomorphism_to_element(MUV, psi)
          push!(img_gens, img_g)
        end

        # Finally, this allows us to assemble the restriction map
        res = hom(MW, MUV, img_gens, OOX(W, UV))
        add_incoming_restriction!(M, W, MUV, res)
      end
    end
  end
  C = Covering(patch_list)
  fill_with_lazy_glueings!(C, scheme(M))
  return C
end

function _trivializing_covering(M::AbsCoherentSheaf, U::AbsSpec)
  X = scheme(M)
  OOX = OO(X)
  MU = M(U)
  MU isa FreeMod && return [U]
  MU::SubQuo
  A = _presentation_matrix(MU)
  if iszero(A) 
    # Trivial shortcut in the recursion. 
    # We nevertheless need to recreate U as a PrincipalOpenSubset of itself 
    # as we are not allowed to alter the values of the sheaf M on U directly.
    V = PrincipalOpenSubset(U, one(OO(U)))
    F = FreeMod(OO(V), ncols(A))
    res = hom(MU, F, gens(F), OOX(U, V))
    add_incoming_restriction(M, U, F, res)
    object_cache(M)[V] = F
    return [V]
  end

  # We do not need to go through all entries of A, but only those 
  # necessary to generate the unit ideal.
  I = ideal(OOX(U), [A[i, j] for i in 1:nrows(A) for j in 1:ncols(A)])
  one(OOX(U)) in I || error("sheaf is not locally trivial")
  # The non-zero coordinates provide us with a list of entries which 
  # are sufficient to do so. This set can not assumed to be minimal, though.
  a = coordinates(one(OOX(U)), I)
  nonzero_entries = [ i for i in 1:ngens(I) if !iszero(a[i])]
  return_patches = AbsSpec[]

  for t in nonzero_entries
    i = div(t-1, ncols(A)) + 1
    j = mod(t-1, ncols(A)) + 1 # The matrix coordinates of the nonzero entry
    # We invert the (i,j)-th entry of A. 
    # Then we can reduce the presentation matrix so that we can throw away one 
    # of the generators of the module. 
    V = PrincipalOpenSubset(U, A[i, j])
    Ares = map_entries(OOX(U, V), A)
    uinv = inv(Ares[i, j])
    multiply_row!(Ares, uinv, i)
    for k in 1:i-1
      #multiply_row!(Ares, u, k)
      add_row!(Ares, -A[k, j], i, k)
    end
    for k in i+1:nrows(Ares)
      #multiply_row!(Ares, u, k)
      add_row!(Ares, -A[k, j], i, k)
    end
    Asub = Ares[[k for k in 1:nrows(Ares) if k != i], [k for k in 1:ncols(Ares) if k !=j]]

    # Assemble the restriction map from the parent node
    if iszero(Asub)
      # End of recursion. 
      # Create a free module and the corresponding restriction morphism.
      F = FreeMod(OO(V), ncols(Asub))
      img_gens = elem_type(F)[]
      for k in 1:j-1
        push!(img_gens, F[k])
      end
      push!(img_gens, 
            -sum([Ares[i, k]*F[(k>j ? k-1 : k)] for k in 1:ncols(Ares) if k!=j], 
                 init=zero(F))
           )
      for k in j+1:ncols(Ares)
        push!(img_gens, F[k-1])
      end
      res = hom(MU, F, img_gens, OOX(U, V))

      # Since we are messing with the internals of the sheaf, we need 
      # to leave everything clean. This includes manual caching.
      add_incoming_restriction!(M, U, F, res)
      object_cache(M)[V] = F
      set_attribute!(F, :_presentation_matrix, Asub)
      push!(return_patches, V)
    else
      # Intermediate recursion step. 
      # Recreate the restriction of the module to the open subset but with one generator 
      # less and construct the restriction map.
      F, amb_res = change_base_ring(OOX(U, V), ambient_free_module(MU))
      v = ambient_representatives_generators(MU)
      M_gens = amb_res.(v)
      rest_gens = [M_gens[k] for k in 1:length(M_gens) if k!=j]
      rels = [amb_res(w) for w in relations(MU)]
      MV = SubQuo(F, rest_gens, rels)
      img_gens = elem_type(F)[]
      for k in 1:j-1
        push!(img_gens, M_gens[k])
      end
      push!(img_gens, 
            -sum([Ares[i, k]*M_gens[k] for k in 1:length(M_gens) if k!=j], 
                 init=zero(F))
           )
      for k in j+1:length(M_gens)
        push!(img_gens, M_gens[k])
      end
      res = hom(MU, MV, MV.(img_gens), OOX(U, V))
      add_incoming_restriction!(M, U, MV, res)
      object_cache(M)[V] = MV
      set_attribute!(MV, :_presentation_matrix, Asub)
      return_patches = vcat(return_patches, _trivializing_covering(M, V))
    end
  end
  return return_patches
end

@attr MatrixElem function _presentation_matrix(M::ModuleFP)
  return matrix(map(presentation(M), 1))
end

function fill_with_lazy_glueings!(C::Covering, X::AbsCoveredScheme)
  # TODO: Fill in all the implicit and restricted glueings from X.
  # We assume that all patches in C are in a tree structure eventually 
  # leading to affine charts of X.
  return C
end

function _common_refinement_list(C::Covering, D::Covering)
  patch_list = AbsSpec[]
  for U in patches(C)
    for V in patches(D)
      success, W = _have_common_ancestor(U, V)
      if success
        #TODO: Model the intersection of both U and V with an appropriate 
        #place in the tree structure.
        incU = _flatten_open_subscheme(U, W)
        incV = _flatten_open_subscheme(V, W)
        UV = intersect(codomain(incU), codomain(incV))
        push!(patch_list, UV)
      end
    end
  end
  return patch_list
end
