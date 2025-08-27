### Lookup and production pattern for sheaves of modules
#
# When asked to produce a module on an open affine U, the functions
# below lead to the following behavior.
#
#                     U₁                    an `affine_chart` of `X`
#           _________/|\______________      (`patches` of `default_covering(X)`)
#          /          |               \
#          V₁         V₂______        |     two `PrincipalOpenSubset`s of U₁
#         /|\        /|\      \       |     covering the latter
#        / | \      / | \     |       |
#       W₁ W₂ W₃   A₁ A₂ A₃   C₁      D₁    `PrincipalOpenSubset`s of, respectively,
#                        |                  V₁, V₂, and U₁
#                        |
#                        E
#
#  Figure 1: A sample tree in one `affine_chart` of a `CoveredScheme`
#
#  Suppose the patches of the `default_covering` of ℱ (i.e. now the local variable
#  `default_cov`) are V₁ together with A₁, A₂, and A₃.
#
# 1) Look up whether U is in the list of patches of `default_cov`, e.g. U = V₁
#    or U = A₂. If yes, return the value of the dictionary.
#
# 2) See whether U hangs below some V in the ancestry tree with V
#    in `default_cov`; e.g. U = E or U = W₃. If yes, restrict from V.
#
# 3) Otherwise, U is covered by patches {Aᵢ}, i ∈ I, from `default_cov` and
#    contained in one affine chart U₁ of X
#      U₁ ⊃ U ⊃ Aᵢ.
#    We distinguish two sub-cases.
#
#    3.1) U appears in the refinement tree T of W whose leafs consist
#    entirely of the Aᵢs and the latter cover U. For instance, this
#    could be the case for U = V₂ with the Aᵢ covering it.
#    Let T' be the subtree starting from U and let {Aⱼ}, j ∈ J be the
#    leafs of this subtree. We then recursively build up the modules on
#    the nodes in T'. Note that this requires some further obvious
#    covering properties on the subtree T'.
#
#    3.2) U does not appear in the refinement tree T above; e.g.
#    when U = C₁ or U = D₁. Then we first have to build the module on
#    U₁ and restrict from there.
#
# The point 3) is not implemented, yet. Instead, the current code
# requires to take refinements hanging under nodes in `default_cov`.

### Production of the modules on open sets; to be cached
function produce_object(
    F::SheafOfModules,
    U::AbsAffineScheme
  )
  MD = F.ID
  MG = F.MG
  X = scheme(F)
  # Since the other cases are caught by the methods below,
  # U can only be an affine_chart of X.
  #
  # See whether we have anything cached for U
  haskey(MD, U) && return MD[U]

  # If not, we are in case 3) above.
  error("production of modules not implemented in this case")
end

function produce_object(
    F::SheafOfModules,
    U::PrincipalOpenSubset
  )
  MD = F.ID
  MG = F.MG
  X = scheme(F)
  # Check whether we can directly produce the module
  haskey(MD, U) && return MD[U]

  # If not: See whether we are below a prescribed module in the
  # refinement tree
  if has_ancestor(x->haskey(MD, x), U)
    V = ambient_scheme(U)
    MV = F(V)
    rho = OO(X)(V, U)
    MU, phi = change_base_ring(rho, MV)
    add_incoming_restriction!(F, V, MU, phi)
    return MU
  end
  # We are in case 3) above
  error("production of modules not implemented in this case")
end

function produce_object(
    F::SheafOfModules,
    U::SimplifiedAffineScheme
  )
  MD = F.ID
  MG = F.MG
  X = scheme(F)
  # Check whether we can directly produce the module
  haskey(MD, U) && return MD[U]

  # If not: See whether we are below a prescribed module in the
  # refinement tree
  if has_ancestor(x->haskey(MD, x), U)
    V = original(U)
    MV = F(V)
    rho = OO(X)(V, U)
    MU, phi = change_base_ring(rho, MV)
    add_incoming_restriction!(F, V, MU, phi)
    return MU
  end
  # We are in case 3) above
  error("production of modules not implemented in this case")
end

### Production of the restriction maps; to be cached
function produce_restriction_map(F::SheafOfModules, V::AbsAffineScheme, U::AbsAffineScheme)
  MD = F.ID
  MG = F.MG
  X = scheme(F)
  # Since the other cases are caught by the methods below, both U and V
  # must be `affine_chart`s of X with U contained in V along some gluing.
  if any(W->(W === U), affine_charts(X)) && any(W->(W === V), affine_charts(X))
    MV = F(V)
    MU = F(U)
    A = MG[(V, U)] # The transition matrix
    UU, _ = gluing_domains(default_covering(X)[U, V])
    psi = OO(X)(UU, U) # Needs to exist by the checks of is_open_func, even though
    # in general UU ⊂ U!
    return hom(MV, MU, [sum([psi(A[i, j]) * MU[j] for j in 1:ngens(MU)], init=zero(MU)) for i in 1:ngens(MV)], OO(X)(V, U))
  else
    error("invalid input")
  end
end

function produce_restriction_map(F::SheafOfModules, V::AbsAffineScheme, U::PrincipalOpenSubset)
  MD = F.ID
  MG = F.MG
  X = scheme(F)
  # If V was not an affine_chart of X, some other function would have
  # been triggered.

  # First the easy case: Inheritance from an ancestor in the tree.
  if ambient_scheme(U) === V
    # If the restriction was more complicated than what follows, then
    # it would have been cached earlier and this call would not have happened
    # This is the end of the recursion induced in the next elseif below.
    res = hom(F(V), F(U), gens(F(U)), OO(X)(W, U))
    return res
  elseif has_ancestor(W->(W === V), U)
    W = ambient_scheme(U)
    return compose(F(V, W), F(W, U))
  end

  # Now we know we have a transition across charts
  W = __find_chart(U, default_covering(X))
  A = MG[(V, W)] # The transition matrix
  WW, _ = gluing_domains(default_covering(X)[W, V])
  # From W to U (and hence also from WW to U) the generators of the modules
  # in F might have changed. Thus, we have to expect a non-trivial transition
  # from the top-level down to U. The transition matrix A is only given with
  # respect to the generators of F(W), so we have to map them manually down.
  # The call to F(W, U) will be handled by the above if-clauses.
  return hom(F(V), F(U),
             [sum([OO(X)(WW, U)(A[i, j])*F(W, U)(F(W)[j]) for j in 1:ngens(F(W))], init=zero(F(U)))
              for i in 1:ngens(F(V))],
             OO(X)(V, U)
            )
end

function produce_restriction_map(F::SheafOfModules, V::AbsAffineScheme, U::SimplifiedAffineScheme)
  MD = F.ID
  MG = F.MG
  X = scheme(F)
  # If V was not an affine_chart of X, some other function would have
  # been triggered.
  @assert any(x->x===V, affine_charts(X)) "first argument must be an affine chart"

  # First the easy case: Inheritance from an ancestor in the tree.
  if original(U) === V
    # If the restriction was more complicated than what follows, then
    # it would have been cached earlier and this call would not have happened
    # This is the end of the recursion induced in the next elseif below.
    W = original(U)
    res = hom(F(W), F(U), gens(MU), OO(X)(W, U))
    return res
  elseif has_ancestor(W->(W === V), U)
    W = original(U)
    return compose(F(V, W), F(W, U))
  end

  # Now we know we have a transition across charts
  W = __find_chart(U, default_covering(X))
  A = MG[(V, W)] # The transition matrix
  WW, _ = gluing_domains(default_covering(X)[W, V])
  # From W to U (and hence also from WW to U) the generators of the modules
  # in F might have changed. Thus, we have to expect a non-trivial transition
  # from the top-level down to U. The transition matrix A is only given with
  # respect to the generators of F(W), so we have to map them manually down.
  # The call to F(W, U) will be handled by the above if-clauses.
  return hom(F(V), F(U),
             [sum([OO(X)(WW, U)(A[i, j])*F(W, U)(F(W)[j]) for j in 1:ngens(F(W))], init=zero(F(U)))
              for i in 1:ngens(F(V))],
             OO(X)(V, U)
            )
end
function produce_restriction_map(F::SheafOfModules,
    V::Union{<:PrincipalOpenSubset, <:SimplifiedAffineScheme},
    U::AbsAffineScheme
  )
  MD = F.ID
  MG = F.MG
  X = scheme(F)
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
            sum(OO(X)(V, U)(c[i])*gens(F(U), i) for i in 1:ngens(F(U)))
           )
    end
    return hom(F(V), F(U), img_gens, OO(X)(V, U))
  else
    # U must be properly contained in the gluing domains of the
    # gluing of the affine chart of V with U.
    error("case not implemented")
  end
  # Problem: We can assume that we know how to pass from generators
  # of W = __find_chart(V, default_covering(X)) to those on V, but we do not
  # know the inverse to this. But the transition matrix to U is given
  # with respect to the generators on W.
  error("case not implemented")
end
function produce_restriction_map(F::SheafOfModules, V::PrincipalOpenSubset, U::PrincipalOpenSubset)
  MD = F.ID
  MG = F.MG
  X = scheme(F)
  V === U && return id_hom(F(U))

  if V === ambient_scheme(U)
    return hom(F(V), F(U), gens(F(U)), OO(X)(V, U)) # If this had been more complicated, it would have been cached.
  elseif has_ancestor(W->W===V, U)
    W = ambient_scheme(U)
    return compose(F(V, W), F(W, U))
  end

  # Below follow the more complicated cases.
  success, _ = _have_common_ancestor(U, V)
  if success
    W = __find_chart(U, default_covering(X))
    gens_U = F(W, U).(gens(F(W))) # This will be caught by the preceding clauses
    gens_V = F(W, V).(gens(F(W)))
    sub_V, inc = sub(F(V), gens_V)
    img_gens = elem_type(F(U))[]
    for v in gens(F(V))
      w = preimage(inc, v) # We know that inc is actually an isomorphism
      c = coordinates(w)
      w = sum(OO(X)(V, U)(c[i])*gens_U[i]
              for i in 1:length(gens_U)
             )
      push!(img_gens, w)
    end
    return hom(F(V), F(U), img_gens, OO(X)(V, U))
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
  gens_U = F(WU, U).(gens(F(WU))) # This will be caught by the preceding clauses
  gens_V = F(WV, V).(gens(F(WV)))
  sub_V, inc = sub(F(V), gens_V)
  img_gens = elem_type(F(U))[]
  A = MG[(WV, WU)] # The transition matrix
  WW, _ = gluing_domains(default_covering(X)[WU, WV])
  for v in gens(F(V))
    w = preimage(inc, v) # We know that inc is actually an isomorphism
    c = coordinates(w)
    w = sum(sum(OO(X)(V, U)(c[i])*OO(X)(WW, U)(A[i, j])*gens_U[j]
                for i in 1:length(gens_V))
            for j in 1:length(gens_U)
           )
    push!(img_gens, w)
  end
  return hom(F(V), F(U), img_gens, OO(X)(V, U))
end
function produce_restriction_map(F::SheafOfModules, V::PrincipalOpenSubset, U::SimplifiedAffineScheme)
  MD = F.ID
  MG = F.MG
  X = scheme(F)
  if V === original(U)
    return hom(F(V), F(U), gens(F(U)), OO(X)(V, U)) # If this had been more complicated, it would have been cached.
  elseif has_ancestor(W->W===V, U)
    W = original(U)
    return compose(F(V, W), F(W, U))
  end

  # Below follow the more complicated cases.
  success, _ = _have_common_ancestor(U, V)
  if success
    W = __find_chart(U, default_covering(X))
    gens_U = F(W, U).(gens(F(W))) # This will be caught by the preceding clauses
    gens_V = F(W, V).(gens(F(W)))
    sub_V, inc = sub(F(V), gens_V)
    img_gens = elem_type(F(U))[]
    for v in gens(F(V))
      w = preimage(inc, v) # We know that inc is actually an isomorphism
      c = coordinates(w)
      w = sum(OO(X)(V, U)(c[i])*gens_U[i]
              for i in 1:length(gens_U)
             )
      push!(img_gens, w)
    end
    return hom(F(V), F(U), img_gens, OO(X)(V, U))
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
  gens_U = F(WU, U).(gens(F(WU))) # This will be caught by the preceding clauses
  gens_V = F(WV, V).(gens(F(WV)))
  sub_V, inc = sub(F(V), gens_V)
  img_gens = elem_type(F(U))[]
  A = MG[(WV, WU)] # The transition matrix
  WW, _ = gluing_domains(default_covering(X)[WU, WV])
  for v in gens(F(V))
    w = preimage(inc, v) # We know that inc is actually an isomorphism
    c = coordinates(w)
    w = sum(sum(OO(X)(V, U)(c[i])*OO(X)(WW, U)(A[i, j])*gens_U[j]
                for i in 1:length(gens_V))
            for j in 1:length(gens_U)
           )
    push!(img_gens, w)
  end
  return hom(F(V), F(U), img_gens, OO(X)(V, U))
end
function produce_restriction_map(F::SheafOfModules, V::SimplifiedAffineScheme, U::PrincipalOpenSubset)
  MD = F.ID
  MG = F.MG
  X = scheme(F)
  if V === ambient_scheme(U)
    return hom(F(V), F(U), gens(F(U)), OO(X)(V, U)) # If this had been more complicated, it would have been cached.
  elseif has_ancestor(W->W===V, U)
    W = ambient_scheme(U)
    return compose(F(V, W), F(W, U))
  end

  # Below follow the more complicated cases.
  success, _ = _have_common_ancestor(U, V)
  if success
    W = __find_chart(U, default_covering(X))
    gens_U = F(W, U).(gens(F(W))) # This will be caught by the preceding clauses
    gens_V = F(W, V).(gens(F(W)))
    sub_V, inc = sub(F(V), gens_V)
    img_gens = elem_type(F(U))[]
    for v in gens(F(V))
      w = preimage(inc, v) # We know that inc is actually an isomorphism
      c = coordinates(w)
      w = sum(OO(X)(V, U)(c[i])*gens_U[i]
              for i in 1:length(gens_U)
             )
      push!(img_gens, w)
    end
    return hom(F(V), F(U), img_gens, OO(X)(V, U))
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
  gens_U = F(WU, U).(gens(F(WU))) # This will be caught by the preceding clauses
  gens_V = F(WV, V).(gens(F(WV)))
  sub_V, inc = sub(F(V), gens_V)
  img_gens = elem_type(F(U))[]
  A = MG[(WV, WU)] # The transition matrix
  WW, _ = gluing_domains(default_covering(X)[WU, WV])
  for v in gens(F(V))
    w = preimage(inc, v) # We know that inc is actually an isomorphism
    c = coordinates(w)
    w = sum(sum(OO(X)(V, U)(c[i])*OO(X)(WW, U)(A[i, j])*gens_U[j]
                for i in 1:length(gens_V))
            for j in 1:length(gens_U)
           )
    push!(img_gens, w)
  end
  return hom(F(V), F(U), img_gens, OO(X)(V, U))
end
function produce_restriction_map(F::SheafOfModules, V::SimplifiedAffineScheme, U::SimplifiedAffineScheme)
  MD = F.ID
  MG = F.MG
  X = scheme(F)
  V === U && return id_hom(F(U))

  if V === original(U)
    return hom(F(V), F(U), gens(F(U)), OO(X)(V, U)) # If this had been more complicated, it would have been cached.
  elseif has_ancestor(W->W===V, U)
    W = original(U)
    return compose(F(V, W), F(W, U))
  end

  # Below follow the more complicated cases.
  success, _ = _have_common_ancestor(U, V)
  if success
    W = __find_chart(U, default_covering(X))
    gens_U = F(W, U).(gens(F(W))) # This will be caught by the preceding clauses
    gens_V = F(W, V).(gens(F(W)))
    sub_V, inc = sub(F(V), gens_V)
    img_gens = elem_type(F(U))[]
    for v in gens(F(V))
      w = preimage(inc, v) # We know that inc is actually an isomorphism
      c = coordinates(w)
      w = sum(OO(X)(V, U)(c[i])*gens_U[i]
              for i in 1:length(gens_U)
             )
      push!(img_gens, w)
    end
    return hom(F(V), F(U), img_gens, OO(X)(V, U))
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
  gens_U = F(WU, U).(gens(F(WU))) # This will be caught by the preceding clauses
  gens_V = F(WV, V).(gens(F(WV)))
  sub_V, inc = sub(F(V), gens_V)
  img_gens = elem_type(F(U))[]
  A = MG[(WV, WU)] # The transition matrix
  WW, _ = gluing_domains(default_covering(X)[WU, WV])
  for v in gens(F(V))
    w = preimage(inc, v) # We know that inc is actually an isomorphism
    c = coordinates(w)
    w = sum(sum(OO(X)(V, U)(c[i])*OO(X)(WW, U)(A[i, j])*gens_U[j]
                for i in 1:length(gens_V))
            for j in 1:length(gens_U)
           )
    push!(img_gens, w)
  end
  return hom(F(V), F(U), img_gens, OO(X)(V, U))
end


### production for HomSheaf
function produce_object(FF::HomSheaf, U::AbsAffineScheme)
  F = domain(FF)
  G = codomain(FF)
  return hom(F(U), G(U))[1]
end

function produce_restriction_map(FF::HomSheaf, V::AbsAffineScheme, U::AbsAffineScheme)
  X = scheme(FF)
  MV = FF(V)
  MU = FF(U)
  F = domain(FF)
  G = codomain(FF)
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

  return hom(MV, MU, images, OO(X)(V, U)) # TODO: Set check=false?
end

### production for DirectSumSheaf
function produce_object(FF::DirectSumSheaf, U::AbsAffineScheme)
  result, inc, pr = direct_sum([F(U) for F in summands(FF)]..., task=:both)
  set_attribute!(result, :inclusions, inc) # TODO: Workaround as long as the maps are not cached.
  set_attribute!(result, :projections, pr)
  return result
end

function produce_restriction_map(FF::DirectSumSheaf, V::AbsAffineScheme, U::AbsAffineScheme)
  X = scheme(FF)
  MV = FF(V)
  MU = FF(U)
  inc_V = get_attribute(MV, :inclusions)::Vector
  pr_V = get_attribute(MV, :projections)::Vector
  inc_U = get_attribute(MU, :inclusions)::Vector
  pr_U = get_attribute(MU, :projections)::Vector

  parts = [] # TODO: Can we do better with type annotation?
  for i in 1:length(inc_V)
    push!(parts, hom(MV, MU,
                     inc_U[i].(summands(FF)[i](V, U).(pr_V[i].(gens(MV)))),
                     OO(X)(V, U)
                    ))
  end
  return sum(parts)
end


### Production for PushforwardSheaf
function produce_object(FF::PushforwardSheaf, U::AbsAffineScheme)
  X = scheme(FF)
  inc = morphism(FF)
  M = original_sheaf(FF)
  ident = identification_dict(FF)
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
function produce_object(FF::PushforwardSheaf, U::PrincipalOpenSubset)
  Y = scheme(FF)
  OOY = OO(Y)
  inc = morphism(FF)
  X = domain(inc)
  M = original_sheaf(FF)
  ident = identification_dict(FF)
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


function produce_restriction_map(FF::PushforwardSheaf, V::AbsAffineScheme, U::AbsAffineScheme)
  Y = scheme(FF)
  OOY = OO(Y)
  inc = morphism(FF)
  X = domain(inc)
  M = original_sheaf(FF)
  ident = identification_dict(FF)
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


### production for PullbackSheaf
#
# Since the morphism f might have an underlying CoveringMorphism ϕ with
# a non-trivial refinement of the `default_covering` of X as a domain,
# we can not expect to easily produce f^*(M) on the `affine_charts` of X.
# Instead, we can produce it on affine opens U ⊂ X which are hanging
# below the patches in `domain(ϕ)`.
#
# For everything else, we proceed as in case 3) of the general
# SheafOfModules, see above.
#
# Again, this case is not implemented for the time being.

function produce_object(FF::PullbackSheaf, U::AbsAffineScheme)
  f = morphism(FF)
  M = original_sheaf(FF)
  OOY = sheaf_of_rings(M)
  X = scheme(FF)
  OOX = OO(X)
  fcov = covering_morphism(f)::CoveringMorphism
  CX = domain(fcov)::Covering
  CY = codomain(fcov)::Covering
  # See whether U is a patch of the domain covering and pull back directly
  if haskey(morphisms(fcov), U)
    floc = morphisms(fcov)[U]
    MU, map = change_base_ring(pullback(floc), M(codomain(floc)))
    pullbacks_on_patches(FF)[U] = map
    return MU
  end

  # We are in case 3).
  error("case not implemented")
end

function produce_object(FF::PullbackSheaf,
    U::Union{<:PrincipalOpenSubset, <:SimplifiedAffineScheme}
  )
  f = morphism(FF)
  M = original_sheaf(FF)
  OOY = sheaf_of_rings(M)
  X = scheme(FF)
  OOX = OO(X)
  fcov = covering_morphism(f)::CoveringMorphism
  CX = domain(fcov)::Covering
  CY = codomain(fcov)::Covering
  # See whether U is a patch of the domain covering and pull back directly
  if haskey(morphisms(fcov), U)
    floc = morphisms(fcov)[U]
    MU, map = change_base_ring(pullback(floc), M(codomain(floc)))
    pullbacks[U] = map
    return MU
  end

  # If not, check whether we are hanging below such a patch in the
  # refinement tree.
  if has_ancestor(y->any(x->(x===y), patches(domain(fcov))), U)
    V = __find_chart(U, domain(fcov))
    MU, res = change_base_ring(OOX(V, U), FF(V))
    add_incoming_restriction!(FF, V, MU, res)
    return MU
  end

  # We are in case 3)
  error("case not implemented")
end

### Restriction for pulled back sheaves of modules
#
# For U ⊂ V ⊂ X, f : X → Y, M on Y and F = f^*M we do the following
# to compute the restriction morphism F(V) → F(U).
# Let ϕ be the `covering_morphism` behind f.
#
#             f : X   →    Y
#
#                 ∪        ∪
#                    ϕ[V]
#    f*ℱ (V)      V   →    V' ↦ ℱ (V')
#
#      ↓ f*ρ      ∪        ∪      ↓ ρ
#                    ϕ[U]
#    f*ℱ (U)      U   →    U' ↦ ℱ (U')
#
# 1) If both U and V are in the `Covering` `domain(ϕ)` induce the restriction
#    f*ρ from ρ on Y.
# 2) If V is in `domain(ϕ)` and U is a node hanging below V in
#    the refinement tree, restrict from V.
# 3) If V is in `domain(ϕ)` and U is a subset of V, restrict as usual.
# 4) If V is a node hanging below some patch in `domain(ϕ)` and
#    U is a subset, restrict as usual.

function restriction_func(F::PullbackSheaf, V::AbsAffineScheme, U::AbsAffineScheme)
  f = morphism(F)
  M = original_sheaf(F)
  OOY = sheaf_of_rings(M)
  OOX = OO(X)
  fcov = covering_morphism(f)::CoveringMorphism
  CX = domain(fcov)::Covering
  CY = codomain(fcov)::Covering
  if haskey(morphisms(fcov), V)
    if haskey(morphisms(fcov), U)
      # case 1)
      f_V = morphisms(fcov)[V]
      f_U = morphisms(fcov)[U]
      MYV = M(codomain(f_V))
      MYU = M(codomain(f_U))
      res_Y = M(codomain(f_V), codomain(f_U))
      result = hom(F(V), F(U),
                   (pullbacks[U]).(res_Y.(gens(MYV))),
                   OOX(V, U))
      return result
    end

    # case 2)
    error("restriction map should have been cached by production")
  end
  error("case not implemented")
end

