export morphism_type, morphisms
export refinements


########################################################################
# Methods for Covering                                                 #
########################################################################


### essential getters


#function add_affine_refinement!(
#    C::Covering, U::SpecOpen; 
#    a::Vector{RingElemType}=as_vector(coordinates(one(OO(ambient_scheme(U))),
#                                                  ideal(OO(ambient_scheme(U)),
#                                                        OO(ambient_scheme(U)).(gens(U)))),
#                                      length(gens(U))),
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

function intersect_in_covering(U::AbsSpec, V::AbsSpec, C::Covering)
  U in C || error("first patch not found in covering")
  V in C || error("second patch not found in covering")
  (i, j, k) = indexin(U, C)
  (l, m, n) = indexin(V, C)
  if i == l # U and V are affine opens of the same patch X in C
    error("affine refinements of coverings not implemented at the moment")
#    X = C[i]
#    if X === U === V
#      iso = identity_map(X)
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
#    isoW = identity_map(W)
#    incWtoU = inclusion_morphism(W, U, check=false)
#    incWtoV = inclusion_morphism(W, V, check=false)
#    return isoW, isoW, incWtoU, incWtoV
  else
    G = C[i, l]
    (f, g) = glueing_morphisms(G)
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
#glueing_type(C::Covering{SpecType, GlueingType, SpecOpenType}) where {SpecType<:Spec, GlueingType<:Glueing, SpecOpenType<:SpecOpen} = GlueingType
#affine_patch_type(::Type{Covering{SpecType, GlueingType, SpecOpenType, RingElemType}}) where {SpecType<:Spec, GlueingType<:Glueing, SpecOpenType<:SpecOpen, RingElemType<:RingElem} = SpecType
#glueing_type(::Type{Covering{SpecType, GlueingType, SpecOpenType}}) where {SpecType<:Spec, GlueingType<:Glueing, SpecOpenType<:SpecOpen} = GlueingType
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

@attr function standard_covering(X::ProjectiveScheme{CRT}) where {CRT<:AbstractAlgebra.Ring}
  CX = affine_cone(X)
  kk = base_ring(X)
  S = ambient_coordinate_ring(X)
  r = fiber_dimension(X)
  U = Vector{AbsSpec}()
  # TODO: Check that all weights are equal to one. Otherwise the routine is not implemented.
  s = symbols(S)
  for i in 0:r
    R, x = PolynomialRing(kk, [Symbol("("*String(s[k+1])*"//"*String(s[i+1])*")") for k in 0:r if k != i])
    phi = hom(S, R, vcat(gens(R)[1:i], [one(R)], gens(R)[i+1:r]))
    I = ideal(R, phi.(gens(defining_ideal(X))))
    push!(U, Spec(quo(R, I)[1]))
  end
  result = Covering(U)
  for i in 1:r
    for j in i+1:r+1
      x = gens(base_ring(OO(U[i])))
      y = gens(base_ring(OO(U[j])))
      Ui = PrincipalOpenSubset(U[i], OO(U[i])(x[j-1]))
      Uj = PrincipalOpenSubset(U[j], OO(U[j])(y[i]))
      f = SpecMor(Ui, Uj,
                      vcat([x[k]//x[j-1] for k in 1:i-1],
                           [1//x[j-1]],
                           [x[k-1]//x[j-1] for k in i+1:j-1],
                           [x[k]//x[j-1] for k in j:r],
                           x[r+1:end]),
                      check=false
                     )
      g = SpecMor(Uj, Ui,
                      vcat([y[k]//y[i] for k in 1:i-1],
                           [y[k+1]//y[i] for k in i:j-2],
                           [1//y[i]],
                           [y[k]//y[i] for k in j:r],
                           y[r+1:end]),
                      check=false
                     )
      add_glueing!(result, SimpleGlueing(U[i], U[j], f, g, check=false))
    end
  end
  return result
end

@attr function standard_covering(X::ProjectiveScheme{CRT}) where {CRT<:Union{<:MPolyQuoLocalizedRing, <:MPolyLocalizedRing, <:MPolyRing, <:MPolyQuo}}
  CX = affine_cone(X)
  Y = base_scheme(X)
  R = ambient_coordinate_ring(Y)
  kk = coefficient_ring(R)
  S = ambient_coordinate_ring(X)
  r = fiber_dimension(X)
  U = Vector{AbsSpec}()
  pU = IdDict{AbsSpec, AbsSpecMor}()
  # TODO: Check that all weights are equal to one. Otherwise the routine is not implemented.
  s = symbols(S)
  # for each homogeneous variable, set up the chart 
  for i in 0:r
    R_fiber, x = PolynomialRing(kk, [Symbol("("*String(s[k+1])*"//"*String(s[i+1])*")") for k in 0:r if k != i])
    F = Spec(R_fiber)
    ambient_space, pF, pY = product(F, Y)
    fiber_vars = pullback(pF).(gens(R_fiber))
    mapped_polys = [map_coefficients(pullback(pY), f) for f in gens(defining_ideal(X))]
    patch = subscheme(ambient_space, [evaluate(f, vcat(fiber_vars[1:i], [one(OO(ambient_space))], fiber_vars[i+1:end])) for f in mapped_polys])
    push!(U, patch)
    pU[patch] = restrict(pY, patch, Y, check=false)
  end
  result = Covering(U)
  for i in 1:r
    for j in i+1:r+1
      x = gens(base_ring(OO(U[i])))
      y = gens(base_ring(OO(U[j])))
      f = SpecOpenMor(U[i], x[j-1], 
                      U[j], y[i],
                      vcat([x[k]//x[j-1] for k in 1:i-1],
                           [1//x[j-1]],
                           [x[k-1]//x[j-1] for k in i+1:j-1],
                           [x[k]//x[j-1] for k in j:r],
                           x[r+1:end]),
                      check=false
                     )
      g = SpecOpenMor(U[j], y[i],
                      U[i], x[j-1],
                      vcat([y[k]//y[i] for k in 1:i-1],
                           [y[k+1]//y[i] for k in i:j-2],
                           [1//y[i]],
                           [y[k]//y[i] for k in j:r],
                           y[r+1:end]),
                      check=false
                     )
      add_glueing!(result, Glueing(U[i], U[j], f, g, check=false))
    end
  end
  covered_projection = CoveringMorphism(result, Covering(Y), pU)
  set_attribute!(X, :covered_projection_to_base, covered_projection)
  return result
end

########################################################################
# Methods for CoveringMorphism                                         #
########################################################################

#covering_type(C::CoveringMorphism{R, S, T}) where {R, S, T} = S
#covering_type(::Type{CoveringMorphism{R, S, T}}) where {R, S, T} = S
#affine_patch_type(C::CoveringMorphism{R, S, T}) where {R, S, T} = R
#affine_patch_type(::Type{CoveringMorphism{R, S, T}}) where {R, S, T} = R

#morphism_type(C::Covering{SpecType, GlueingType, SpecOpenType}) where {SpecType<:Spec, GlueingType<:Glueing, SpecOpenType<:SpecOpen} = CoveringMorphism{SpecType, Covering{SpecType, GlueingType, SpecOpenType}, morphism_type(SpecType, SpecType)}
#morphism_type(::Type{Covering{SpecType, GlueingType, SpecOpenType}}) where {SpecType<:Spec, GlueingType<:Glueing, SpecOpenType<:SpecOpen} = CoveringMorphism{SpecType, Covering{SpecType, GlueingType, SpecOpenType}, morphism_type(SpecType, SpecType)}



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
#covered_scheme_type(::Type{T}) where {T<:Spec} = CoveredScheme{covering_type(T), morphism_type(covering_type(T))}
#covered_scheme_type(X::Spec) = covered_scheme_type(typeof(X))
#
#covered_scheme_type(::Type{T}) where {T<:ProjectiveScheme} = covered_scheme_type(affine_patch_type(P))
#covered_scheme_type(P::ProjectiveScheme) = covered_scheme_type(typeof(P))

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
#@Markdown.doc """
#    common_refinement(X::CoveredScheme, C1::T, C2::T) where {T<:Covering}
#
#Given two coverings of ``X``, return a triple `(C_new, f, g)` consisting 
#of a common refinement `C_new` of `C1` and `C2` and the refinement morphisms 
#`f : C_new → C1` and `g : C_new → C2`.
#"""
#function common_refinement(X::CoveredScheme, C1::T, C2::T) where {T<:Covering}
#  # shortcut for the trivial cases
#  C1 == C2 && return (C1, identity_map(C1), identity_map(C1))
#
#  # find the minimal common root using the refinement graph
#  r, p1, p2 = find_common_root(refinement_graph(X), X[C1], X[C2])
#
#  # if one covering sits strictly on top of the other, take the shortcut
#  if length(p1) == 0
#    return (C2, identity_map(C1), _compose_along_path(X, p2))
#  end
#  if length(p2) == 0
#    return (C1, _compose_along_path(X, p1), identity_map(C2))
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
#      if issubset(U, V) 
#        inc1[U] = identity_map(U)
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
#  # cook up the glueings for the new patches from those in the common root.
#  new_glueings = IdDict{Tuple{affine_patch_type(X), affine_patch_type(X)}, glueing_type(affine_patch_type(X))}()
#  for (W1, W2) in keys(glueings(C0))
#    U_patches = [U for U in new_patches if codomain(inc0[U]) === W1]
#    V_patches = [V for V in new_patches if codomain(inc0[V]) === W2]
#    for U in U_patches
#      for V in V_patches
#        new_glueings[(U, V)] = restrict(C0[W1, W2], U, V)
#      end
#    end
#  end
#
#  C_new = Covering(new_patches, new_glueings)
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
