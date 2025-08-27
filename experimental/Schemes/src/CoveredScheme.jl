export morphism_type

########################################################################
# Methods for Covering                                                 #
########################################################################

### essential getters

#function add_affine_refinement!(
#    C::Covering, U::AffineSchemeOpenSubscheme;
#    a::Vector{RingElemType}=as_vector(coordinates(one(OO(ambient_scheme(U))),
#                                                  ideal(OO(ambient_scheme(U)),
#                                                        OO(ambient_scheme(U)).(gens(U)))),
#                                      ngens(U)),
#    check::Bool=true
#  ) where {RingElemType<:RingElem}
#  X = ambient_scheme(U)
#  @show typeof(OO(X))
#  X in patches(C) || error("ambient scheme not found in the basic patches of the covering")
#  if check
#    isone(OO(X)(sum([c*g for (c, g) in zip(a, gens(U))]))) || error("patches of $U do not cover $X")
#  end
#  if !haskey(affine_refinements(C), X)
#    affine_refinements(C)[X] = [(U, a)]
#  else
#    push!(affine_refinements(C)[X], (U, a))
#  end
#  return C
#end

# functions for handling sets in coverings

function intersect_in_covering(U::AbsAffineScheme, V::AbsAffineScheme, C::Covering)
  U in C || error("first patch not found in covering")
  V in C || error("second patch not found in covering")
  (i, j, k) = indexin(U, C)
  (l, m, n) = indexin(V, C)
  if i == l # U and V are affine opens of the same patch X in C
    error("affine refinements of coverings not implemented at the moment")
#    X = C[i]
#    if X === U === V
#      iso = id_hom(X)
#      return iso, iso, iso, iso
#    end
#    f = one(ambient_coordinate_ring(X)) # For the case where U is already a basic patch
#    if j != 0 # In this case, U appears in some refinement
#      f = gens(affine_refinements(C)[X][j][1])[k]
#    end
#    g = one(ambient_coordinate_ring(X))
#    if m != 0
#      g = gens(affine_refinements(C)[X][m][1])[n]
#    end
#    W = PrincipalOpenSubset(X, [f*g])
#    isoW = id_hom(W)
#    incWtoU = inclusion_morphism(W, U, check=false)
#    incWtoV = inclusion_morphism(W, V, check=false)
#    return isoW, isoW, incWtoU, incWtoV
  else
    G = C[i, l]
    (f, g) = gluing_morphisms(G)
    preimV = preimage(f, V)
    preimU = preimage(g, U)
    WU = intersect(preimV, U)
    WV = intersect(preimU, V)
    isoWUtoWV = restrict(f, WU, WV, check=false)
    isoWVtoWU = restrict(g, WV, WU, check=false)
    incWUtoU = inclusion_morphism(WU, U, check=false)
    incWVtoV = inclusion_morphism(WV, V, check=false)
    return isoWUtoWV, isoWVtoWU, incWUtoU, incWVtoV
  end
end

#affine_patch_type(C::Covering) = affine_patch_type(typeof(C))
#gluing_type(C::Covering{AffineSchemeType, GluingType, AffineSchemeOpenSubschemeType}) where {AffineSchemeType<:AffineScheme, GluingType<:Gluing, AffineSchemeOpenSubschemeType<:AffineSchemeOpenSubscheme} = GluingType
#affine_patch_type(::Type{Covering{AffineSchemeType, GluingType, AffineSchemeOpenSubschemeType, RingElemType}}) where {AffineSchemeType<:AffineScheme, GluingType<:Gluing, AffineSchemeOpenSubschemeType<:AffineSchemeOpenSubscheme, RingElemType<:RingElem} = AffineSchemeType
#gluing_type(::Type{Covering{AffineSchemeType, GluingType, AffineSchemeOpenSubschemeType}}) where {AffineSchemeType<:AffineScheme, GluingType<:Gluing, AffineSchemeOpenSubschemeType<:AffineSchemeOpenSubscheme} = GluingType
#open_subset_type(::Type{Covering{R, S, T}}) where {R, S, T} = T
#open_subset_type(C::Covering) = open_subset_type(typeof(C))

# TODO: For some reason, the `indexin` method won't work. In the long
# run, one should probably find out why and fix it.


affine_refinements(C::Covering) = C.affine_refinements

### type getters
# Not required at the moment; should eventually be deleted.
#base_morphism_type(::Type{T}) where {DT, CT, BMT, T<:CoveringMorphism{DT, CT, BMT}} = BMT
#base_morphism_type(C::Covering) = base_morphism_type(typeof(C))

#domain_type(::Type{T}) where {DT, CT, BMT, T<:CoveringMorphism{DT, CT, BMT}} = DT
#domain_type(C::Covering) = domain_type(typeof(C))

#codomain_type(::Type{T}) where {DT, CT, BMT, T<:CoveringMorphism{DT, CT, BMT}} = CT
#codomain_type(C::Covering) = codomain_type(typeof(C))

########################################################################
# Constructors for standard schemes (Projective space, etc.)           #
########################################################################

@doc raw"""
    _generate_affine_charts(X::Scheme) -> Dict{Int, AbsAffineScheme}

Helper to generate the affine charts of projective space for `standard_covering`.
This should be overwritten if you want your charts to be of a type different from `AffineScheme`,
for instance `AffinePlaneCurve`.
"""
_generate_affine_charts(X::Scheme)

# The case of a non-trivial homogeneous modulus
function _generate_affine_charts(X::AbsProjectiveScheme{<:Ring, <:MPolyQuoRing})
  chart_dict = Dict{Int, AffineScheme}()
  kk = base_ring(X)
  S = ambient_coordinate_ring(X)
  r = relative_ambient_dimension(X)
  s = symbols(S)
  for i in 0:r
    R, gens_R = polynomial_ring(kk, [Symbol("(",s[k+1],"//",s[i+1],")") for k in 0:r if k != i])
    phi = hom(S, R, vcat(gens_R[1:i], [one(R)], gens_R[i+1:r]); check=false)
    I = ideal(R, phi.(gens(defining_ideal(X))))
    if !isone(I) # return the non-empty charts only
      chart_dict[i+1] = spec(quo(R, I)[1])
    end
  end
  return chart_dict
end

# The case of a trivial homogeneous modulus
function _generate_affine_charts(X::AbsProjectiveScheme{<:Ring, <:MPolyDecRing})
  chart_dict = Dict{Int, AffineScheme}()
  kk = base_ring(X)
  S = ambient_coordinate_ring(X)
  r = relative_ambient_dimension(X)
  s = symbols(S)
  for i in 0:r
    R, x = polynomial_ring(kk, [Symbol("(",s[k+1],"//",s[i+1],")") for k in 0:r if k != i])
    chart_dict[i+1] = spec(R)
  end
  return chart_dict
end

# The case of projective plane curves which wraps the result in the corresponding type
function _generate_affine_charts(X::ProjectivePlaneCurve)
  prelim_dict = _generate_affine_charts(underlying_scheme(X))
  # TODO: Do we want to do this manually and also set the defining_equation?
  return Dict{Int, AffinePlaneCurve}(i => AffinePlaneCurve(U; check=false) for (i, U) in prelim_dict)
end

# The two cases for non-trivial base rings
function _generate_affine_charts(X::AbsProjectiveScheme{<:CRT, <:MPolyQuoRing}) where {CRT<:Union{<:MPolyQuoLocRing, <:MPolyLocRing, <:MPolyRing, <:MPolyQuoRing}}

  chart_dict = Dict{Int, AbsAffineScheme}()
  Y = base_scheme(X)
  R = ambient_coordinate_ring(Y)
  kk = coefficient_ring(R)
  S = ambient_coordinate_ring(X)
  s = symbols(S)
  pU = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()

  r = relative_ambient_dimension(X)
  for i in 0:r
    R_fiber, gens_R_fiber = polynomial_ring(kk, [Symbol("(",s[k+1],"//",s[i+1],")") for k in 0:r if k != i])
    F = spec(R_fiber)
    ambient_space, pF, pY = product(F, Y)
    fiber_vars = pullback(pF).(gens_R_fiber)
    mapped_polys = [map_coefficients(pullback(pY), f) for f in gens(defining_ideal(X))]
    patch = subscheme(ambient_space, elem_type(OO(ambient_space))[evaluate(f, vcat(fiber_vars[1:i], [one(OO(ambient_space))], fiber_vars[i+1:end])) for f in mapped_polys])
    chart_dict[i+1] = patch
    pU[patch] = restrict(pY, patch, Y, check=false)
  end

  return chart_dict, pU
end

function _generate_affine_charts(X::AbsProjectiveScheme{<:CRT, <:MPolyDecRing}) where {CRT<:Union{<:MPolyQuoLocRing, <:MPolyLocRing, <:MPolyRing, <:MPolyQuoRing}}

  chart_dict = Dict{Int, AbsAffineScheme}()
  Y = base_scheme(X)
  R = ambient_coordinate_ring(Y)
  kk = coefficient_ring(R)
  S = ambient_coordinate_ring(X)
  s = symbols(S)
  pU = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()

  r = relative_ambient_dimension(X)
  for i in 0:r
    R_fiber, gens_R_fiber = polynomial_ring(kk, [Symbol("(",s[k+1],"//",s[i+1],")") for k in 0:r if k != i])
    F = spec(R_fiber)
    ambient_space, pF, pY = product(F, Y)
    fiber_vars = pullback(pF).(gens_R_fiber)
    chart_dict[i+1] = ambient_space
    pU[ambient_space] = pY
  end

  return chart_dict, pU
end

# coefficient ring a field, the integers, or similar
@attr Covering function standard_covering(X::AbsProjectiveScheme{<:Ring, <:Union{<:MPolyQuoRing, <:MPolyDecRing}})
  kk = base_ring(X)
  S = ambient_coordinate_ring(X)
  r = relative_ambient_dimension(X)
  # TODO: Check that all weights are equal to one. Otherwise the routine is not implemented.
  s = symbols(S)
  chart_dict = _generate_affine_charts(X)
  iszero(length(chart_dict)) && return empty_covering(base_ring(X))
  @assert all(k->!isempty(chart_dict[k]), keys(chart_dict)) "empty chart created"
  decomp_info = IdDict{AbsAffineScheme, Vector{RingElem}}()
  for (i, U) in chart_dict
    decomp_info[U] = gens(OO(U))[1:i-1]
    _dehomogenization_cache(X)[U] = _dehomogenization_map(X, U, i)
    _homogenization_cache(X)[U] = _homogenization_map(X, U, i)
  end
  result = Covering([chart_dict[i] for i in sort(collect(keys(chart_dict)))])
  set_decomposition_info!(result, decomp_info)
  for (i, U) in chart_dict
    for (j, V) in chart_dict
      j <= i && continue # TODO: Is there a better way to go through this?
      x = gens(ambient_coordinate_ring(U))
      y = gens(ambient_coordinate_ring(V))
      Ui = PrincipalOpenSubset(U, OO(U)(x[j-1]))
      Uj = PrincipalOpenSubset(V, OO(V)(y[i]))
      imgs_f = vcat([x[k]//x[j-1] for k in 1:i-1],
                  [1//x[j-1]],
                  [x[k-1]//x[j-1] for k in i+1:j-1],
                  [x[k]//x[j-1] for k in j:r],
                  x[r+1:end])
      f = morphism(Ui, Uj, [OO(Ui)(a, check=false) for a in imgs_f], check=false)
      imgs_g = vcat([y[k]//y[i] for k in 1:i-1],
                    [y[k+1]//y[i] for k in i:j-2],
                    [1//y[i]],
                    [y[k]//y[i] for k in j:r],
                    y[r+1:end])
      g = morphism(Uj, Ui, [OO(Uj)(b, check=false) for b in imgs_g], check=false)
      add_gluing!(result, SimpleGluing(U, V, f, g, check=false))
    end
  end
  return result
end

# coefficient ring an MPolyAnyRing
@attr Any function standard_covering(X::AbsProjectiveScheme{CRT, <:Union{<:MPolyDecRing, <:MPolyQuoRing}}) where {CRT<:Union{<:MPolyQuoLocRing, <:MPolyLocRing, <:MPolyRing, <:MPolyQuoRing}}
  Y = base_scheme(X)
  R = ambient_coordinate_ring(Y)
  kk = coefficient_ring(R)
  S = ambient_coordinate_ring(X)
  r = relative_ambient_dimension(X)
  U = Vector{AbsAffineScheme}()

  # The case of ℙ⁰-bundles appears frequently in blowups when the
  # ideal sheaf is trivial on some affine open part.
  if r == 0
    result = Covering(Y)
    set_decomposition_info!(result, Y, elem_type(OO(Y))[])
    pU = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
    pU[Y] = id_hom(Y)
    covered_projection = CoveringMorphism(result, result, pU, check=false)
    set_attribute!(X, :covering_projection_to_base, covered_projection)
    _dehomogenization_cache(X)[Y] = _dehomogenization_map(X, Y, 1)
    _homogenization_cache(X)[Y] = _homogenization_map(X, Y, 1)
    return result
  end

  # TODO: Check that all weights are equal to one. Otherwise the routine is not implemented.
  s = symbols(S)
  # for each homogeneous variable, set up the chart
  chart_dict, projection_dict = _generate_affine_charts(X)
  isempty(chart_dict) && return empty_covering(base_ring(Y))

  decomp_info = IdDict{AbsAffineScheme, Vector{RingElem}}()
  for (i, U) in chart_dict
    decomp_info[U] = gens(OO(U))[1:i-1]
    _dehomogenization_cache(X)[U] = _dehomogenization_map(X, U, i)
    _homogenization_cache(X)[U] = _homogenization_map(X, U, i)
  end
  result = Covering([chart_dict[i] for i in sort(collect(keys(chart_dict)))])
  set_decomposition_info!(result, decomp_info)

  for (i, U) in chart_dict
    for (j, V) in chart_dict
      j <= i && continue

      x = gens(ambient_coordinate_ring(U))
      y = gens(ambient_coordinate_ring(V))
      Ui = PrincipalOpenSubset(U, OO(U)(x[j-1]))
      Uj = PrincipalOpenSubset(V, OO(V)(y[i]))
      imgs_f = vcat([x[k]//x[j-1] for k in 1:i-1],
                  [1//x[j-1]],
                  [x[k-1]//x[j-1] for k in i+1:j-1],
                  [x[k]//x[j-1] for k in j:r],
                  x[r+1:end])
      f = morphism(Ui, Uj, [OO(Ui)(a, check=false) for a in imgs_f], check=false)
      imgs_g = vcat([y[k]//y[i] for k in 1:i-1],
                    [y[k+1]//y[i] for k in i:j-2],
                    [1//y[i]],
                    [y[k]//y[i] for k in j:r],
                    y[r+1:end])
      g = morphism(Uj, Ui, [OO(Uj)(b, check=false) for b in imgs_g], check=false)
      add_gluing!(result, SimpleGluing(U, V, f, g, check=false))
    end
  end
  covered_projection = CoveringMorphism(result, Covering(Y), projection_dict, check=false)
  set_attribute!(X, :covering_projection_to_base, covered_projection)
  return result
end


########################################################################
# Methods for CoveringMorphism                                         #
########################################################################

#covering_type(C::CoveringMorphism{R, S, T}) where {R, S, T} = S
#covering_type(::Type{CoveringMorphism{R, S, T}}) where {R, S, T} = S
#affine_patch_type(C::CoveringMorphism{R, S, T}) where {R, S, T} = R
#affine_patch_type(::Type{CoveringMorphism{R, S, T}}) where {R, S, T} = R

#morphism_type(C::Covering{AffineSchemeType, GluingType, AffineSchemeOpenSubschemeType}) where {AffineSchemeType<:AffineScheme, GluingType<:Gluing, AffineSchemeOpenSubschemeType<:AffineSchemeOpenSubscheme} = CoveringMorphism{AffineSchemeType, Covering{AffineSchemeType, GluingType, AffineSchemeOpenSubschemeType}, morphism_type(AffineSchemeType, AffineSchemeType)}
#morphism_type(::Type{Covering{AffineSchemeType, GluingType, AffineSchemeOpenSubschemeType}}) where {AffineSchemeType<:AffineScheme, GluingType<:Gluing, AffineSchemeOpenSubschemeType<:AffineSchemeOpenSubscheme} = CoveringMorphism{AffineSchemeType, Covering{AffineSchemeType, GluingType, AffineSchemeOpenSubschemeType}, morphism_type(AffineSchemeType, AffineSchemeType)}


refinements(X::AbsCoveredScheme) = refinements(underlying_scheme(X))::Dict{<:Tuple{<:Covering, <:Covering}, <:CoveringMorphism}

########################################################################
# Methods for CoveredScheme                                            #
########################################################################
### type getters
#covering_type(X::CoveredScheme{S, T}) where {S, T} = S
#covering_type(::Type{CoveredScheme{S, T}}) where {S, T} = S
#covering_morphism_type(X::CoveredScheme{S, T}) where {S, T} = T
#covering_morphism_type(::Type{CoveredScheme{S, T}}) where {S, T} = T
#affine_patch_type(X::CoveredSchemeType) where {CoveredSchemeType<:CoveredScheme} = affine_patch_type(covering_type(CoveredSchemeType))
#affine_patch_type(::Type{CoveredSchemeType}) where {CoveredSchemeType<:CoveredScheme} = affine_patch_type(covering_type(CoveredSchemeType))

### type constructors
#covered_scheme_type(::Type{T}) where {T<:AffineScheme} = CoveredScheme{covering_type(T), morphism_type(covering_type(T))}
#covered_scheme_type(X::AffineScheme) = covered_scheme_type(typeof(X))
#
#covered_scheme_type(::Type{T}) where {T<:ProjectiveScheme} = covered_scheme_type(affine_patch_type(P))
#covered_scheme_type(P::AbsProjectiveScheme) = covered_scheme_type(typeof(P))

### getter methods
refinements(X::CoveredScheme) = X.refinements

#function set_default_covering!(X::CoveredScheme, C::Covering)
#  C in coverings(X) || error("covering is not listed")
#  X.default_covering = C
#  return X
#end


_compose_along_path(X::CoveredScheme, p::Vector{Int}) = _compose_along_path(X, [X[i] for i in p])

#function _compose_along_path(X::CoveredScheme, p::Vector{CoveringType}) where {CoveringType<:Covering}
#  root = pop!(p)
#  next = pop!(p)
#  mor = X[next, root]
#  while length(p) > 0
#    leaf = pop!(p)
#    mor = compose(X[leaf, next], mor)
#    next = leaf
#  end
#  X[leaf, root] = mor
#  add_edge!(refinement_graph(X), X[leaf], X[root])
#  return mor
#end
#
## TODO: Replace by the polymake routines, once provided!
#function find_common_root(G::Graph{Directed}, i::Int, j::Int)
#  p = [i]
#  Ni = neighbors(G, i)
#  while length(Ni) > 0
#    push!(p, Ni[1])
#    Ni = neighbors(G, Ni[1])
#  end
#  q = [j]
#  Nj = neighbors(G, j)
#  while length(Nj) > 0
#    push!(p, Nj[1])
#    Nj = neighbors(G, Nj[1])
#  end
#  last(p) == last(q) || error("no common root found")
#  return last(p), p, q
#end
#
#@doc raw"""
#    common_refinement(X::CoveredScheme, C1::T, C2::T) where {T<:Covering}
#
#Given two coverings of ``X``, return a triple `(C_new, f, g)` consisting
#of a common refinement `C_new` of `C1` and `C2` and the refinement morphisms
#`f : C_new → C1` and `g : C_new → C2`.
#"""
#function common_refinement(X::CoveredScheme, C1::T, C2::T) where {T<:Covering}
#  # shortcut for the trivial cases
#  C1 == C2 && return (C1, id_hom(C1), id_hom(C1))
#
#  # find the minimal common root using the refinement graph
#  r, p1, p2 = find_common_root(refinement_graph(X), X[C1], X[C2])
#
#  # if one covering sits strictly on top of the other, take the shortcut
#  if length(p1) == 0
#    return (C2, id_hom(C1), _compose_along_path(X, p2))
#  end
#  if length(p2) == 0
#    return (C1, _compose_along_path(X, p1), id_hom(C2))
#  end
#
#  # now we may assume that neither one of the coverings is contained in the other
#  C0 = X[r]
#  f = _compose_along_path(X, p1)
#  g = _compose_along_path(X, p2)
#
#  # prepare for the common refinement
#  new_patches = Vector{affine_patch_type(X)}()
#  inc1 = IdDict{affine_patch_type(X), morphism_type(affine_patch_type(X))}()
#  inc2 = IdDict{affine_patch_type(X), morphism_type(affine_patch_type(X))}()
#  inc0 = IdDict{affine_patch_type(X), morphism_type(affine_patch_type(X))}()
#  for U in patches(C1)
#    W = codomain(f[U])
#    V_candidates = [V for V in patches(C2) if codomain(g[V]) === W]
#
#    # first try to find a patch in C2 which fully includes U
#    patch_found = false
#    while length(V_candidates) > 0
#      V = pop!(V_candidates)
#      if is_subscheme(U, V)
#        inc1[U] = id_hom(U)
#        inc2[U] = inclusion_morphism(U, V)
#        inc0[U] = f[U]
#        push!(new_patches, U)
#        patch_found = true
#        break
#      end
#    end
#    patch_found && break
#
#    # this is the worst case where there is no patch in C2 containing U.
#    for V in V_candidates
#      UV = intersect(U, V)
#      inc1[UV] = inclusion_morphism(UV, U)
#      inc2[UV] = inclusion_morphism(UV, V)
#      inc0[UV] = inclusion_morphism(UV, W)
#      push!(new_patches, UV)
#    end
#  end
#
#  # cook up the gluings for the new patches from those in the common root.
#  new_gluings = IdDict{Tuple{affine_patch_type(X), affine_patch_type(X)}, gluing_type(affine_patch_type(X))}()
#  for (W1, W2) in keys(gluings(C0))
#    U_patches = [U for U in new_patches if codomain(inc0[U]) === W1]
#    V_patches = [V for V in new_patches if codomain(inc0[V]) === W2]
#    for U in U_patches
#      for V in V_patches
#        new_gluings[(U, V)] = restrict(C0[W1, W2], U, V)
#      end
#    end
#  end
#
#  C_new = Covering(new_patches, new_gluings)
#  f = CoveringMorphism(C_new, C1, inc1, check=true) # set to false after debugging
#  g = CoveringMorphism(C_new, C2, inc2, check=true)
#  h = CoveringMorphism(C_new, C0, inc0, check=true)
#  X[C_new, C1] = f
#  X[C_new, C2] = g
#  X[C_new, C0] = h
#  add_edge!(refinement_graph(X), X[C_new], X[C1])
#  add_edge!(refinement_graph(X), X[C_new], X[C2])
#  add_edge!(refinement_graph(X), X[C_new], X[C0])
#  return (C_new, f, g)
#end


### Miscellaneous helper routines
#function as_vector(v::SRow{T}, n::Int) where {T<:RingElem}
#  R = base_ring(v)
#  result = elem_type(R)[zero(R) for i in 1:n]
#  for (i, a) in v
#    result[i] = a
#  end
#  return result
#end



########################################################################
# Bound functionality
########################################################################
function pushforward(f::CompositeCoveredSchemeMorphism, a::VarietyFunctionFieldElem)
  result = a
  for g in maps(f)
    result = pushforward(g, result)
  end
  return result
end

function pullback(f::CompositeCoveredSchemeMorphism, a::VarietyFunctionFieldElem)
  result = a
  for g in reverse(maps(f))
    result = pullback(g, result)
  end
  return result
end

