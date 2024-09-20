#######################################
### Forget toric structure
#######################################

@doc raw"""
    forget_toric_structure(X::AffineNormalToricVariety)

Return a pair `(Y, iso)` where `Y` is a scheme without toric structure,
together with an isomorphism `iso : Y → X`.

# Examples
```jldoctest
julia> C = positive_hull([1 0; 0 1])
Polyhedral cone in ambient dimension 2

julia> antv = affine_normal_toric_variety(C)
Normal toric variety

julia> forget_toric_structure(antv)
(scheme(0), Hom: scheme(0) -> normal toric variety)
```
"""
function forget_toric_structure(X::AffineNormalToricVariety)
  Y = underlying_scheme(X)
  iso = morphism(Y, X, identity_map(OO(X)), check=true)
  iso_inv = morphism(X, Y, identity_map(OO(X)), check=true)
  set_attribute!(iso, :inverse => iso_inv)
  set_attribute!(iso_inv, :inverse => iso)
  return Y, iso
end

@doc raw"""
    forget_toric_structure(X::NormalToricVariety)

Return a pair `(Y, iso)` where `Y` is a scheme without toric structure,
together with an isomorphism `iso : Y → X`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> forget_toric_structure(P2)
(Scheme over QQ covered with 3 patches, Hom: scheme over QQ covered with 3 patches -> normal toric variety)
```
"""
function forget_toric_structure(X::NormalToricVariety)
  # Collect all the isomorphisms forgetting the toric structure
  iso_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
  for U in affine_charts(X)
    iso_dict[U] = forget_toric_structure(U)[2] # store only the isomorphism
  end
  cov = Covering([domain(phi) for (U, phi) in iso_dict])

  # Prepare a dictionary that can be used in the constructor of the covering morphism
  iso_dict_covariant = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
  for (U, phi) in iso_dict
    iso_dict_covariant[domain(phi)] = phi
  end

  # Recreate all the gluings with the new patches
  for U in affine_charts(X), V in affine_charts(X)
    glue = default_covering(X)[U, V]
    new_glue = restrict(glue, inverse(iso_dict[U]), inverse(iso_dict[V]), check=true)
    add_gluing!(cov, new_glue)
  end

  # Prepare the underlying covering morphisms for the identifying isomorphisms
  iso_cov = CoveringMorphism(cov, default_covering(X), iso_dict_covariant)
  inv_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
  for (U, phi) in iso_dict
    inv_dict[U] = inverse(phi)
  end
  iso_cov_inv = CoveringMorphism(default_covering(X), cov, inv_dict)

  # Create the actual scheme without toric structure
  Y = CoveredScheme(cov)

  # Make the identifying isomorphisms
  iso = CoveredSchemeMorphism(Y, X, iso_cov)
  iso_inv = CoveredSchemeMorphism(X, Y, iso_cov_inv)

  set_attribute!(iso, :inverse => iso_inv)
  set_attribute!(iso_inv, :inverse => iso)

  return Y, iso
end


#######################################
### Underlying scheme
#######################################

@doc raw"""
    underlying_scheme(X::AffineNormalToricVariety)

For an affine toric scheme ``X``, this returns
the underlying scheme. In other words, by applying
this method, you obtain a scheme that has forgotten
its toric origin.

# Examples
```jldoctest
julia> C = positive_hull([1 0; 0 1])
Polyhedral cone in ambient dimension 2

julia> antv = affine_normal_toric_variety(C)
Normal toric variety

julia> Oscar.underlying_scheme(antv)
Spectrum
  of quotient
    of multivariate polynomial ring in 2 variables x1, x2
      over rational field
    by ideal (0)
```
"""
@attr AffineScheme{QQField, MPolyQuoRing{QQMPolyRingElem}} underlying_scheme(X::AffineNormalToricVariety) = spec(base_ring(toric_ideal(X)), toric_ideal(X))

###
# Some additional structure to make computation of toric gluings lazy
struct ToricGluingData
  X::NormalToricVariety
  U::AffineNormalToricVariety
  V::AffineNormalToricVariety
# i::Int # The indices of the charts to be glued
# j::Int
end

function _compute_gluing(gd::ToricGluingData)
  X = gd.X
  U = gd.U
  V = gd.V
# i = gd.i
# j = gd.j
# U = affine_charts(X)[i]
# V = affine_charts(X)[j]
  sigma_1 = cone(U)
  sigma_2 = cone(V)
  tau = intersect(sigma_1, sigma_2)
  sigma_1_dual = weight_cone(U)
  sigma_2_dual = weight_cone(V)
  tau_dual = polarize(tau)

  # We do the following. There is a commutative diagram of rings
  #
  #       ℚ [σ₁̌] ↪  ℚ [τ ̌]  ↩  ℚ [σ₂̌]
  #
  # given by localization maps. The cone τ ̌ has lineality L.
  # We need to find a Hilbert basis for both L ∩ σ₁̌ and L ∩ σ₂̌.
  # Then the localization maps are given by inverting the
  # elements of these Hilbert bases. The gluing isomorphisms
  # are then obtained by expressing the generators on the one
  # side in terms of the others.

  # We are using Proposition 1.2.10 in Cox-Little-Schenck here:
  #  "If τ is a face of a polyhedral cone σ and τ* = σ ̌ ∩ τ⟂,
  #   then τ* is a face of σ ̌."

  degs1 = hilbert_basis(U)
  non_local_indices_1 = filter(i->!(vec(-degs1[i,:]) in tau_dual), 1:nrows(degs1))
  degs2 = hilbert_basis(V)
  non_local_indices_2 = filter(i->!(vec(-degs2[i,:]) in tau_dual), 1:nrows(degs2))

  x = gens(OO(U))
  UV = PrincipalOpenSubset(U, [x[i] for i in 1:length(x) if !(i in non_local_indices_1)])
  y = gens(OO(V))
  VU = PrincipalOpenSubset(V, [y[i] for i in 1:length(y) if !(i in non_local_indices_2)])

  y_to_x = _convert_degree_system(degs1, degs2, non_local_indices_1)
  x_to_y = _convert_degree_system(degs2, degs1, non_local_indices_2)

  xx = gens(OO(UV))
  yy = gens(OO(VU))
  f = morphism(UV, VU, [prod((e[i] >= 0 ? u^e[i] : inv(u)^-e[i]) for (i, u) in enumerate(xx); init=one(OO(UV))) for e in y_to_x], check=false)
  g = morphism(VU, UV, [prod((e[i] >= 0 ? v^e[i] : inv(v)^-e[i]) for (i, v) in enumerate(yy); init=one(OO(VU))) for e in x_to_y], check=false)
  set_attribute!(f, :inverse, g)
  set_attribute!(g, :inverse, f)

  result = Gluing(U, V, f, g, check=false)
  return result
end

# Write the elements in `degs2` as linear combinations of `degs1`, allowing only non-negative
# coefficients for the vectors vᵢ of `degs1` with index i ∈ `non_local_indices`.
function _convert_degree_system(degs1::ZZMatrix, degs2::ZZMatrix, non_local_indices_1::Vector{Int})
  result = Vector{ZZMatrix}()
  for i in 1:nrows(degs2)
    C = identity_matrix(ZZ, nrows(degs1))[non_local_indices_1,:]
    S = solve_mixed(transpose(degs1), transpose(degs2[i:i,:]), C; permit_unbounded=true)
    push!(result, S[1:1, :])
  end
  return result
end


@doc raw"""
    underlying_scheme(X::NormalToricVariety)

For a toric covered scheme ``X``, this returns
the underlying scheme. In other words, by applying
this method, you obtain a scheme that has forgotten
its toric origin.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> Oscar.underlying_scheme(P2)
Scheme
  over rational field
with default covering
  described by patches
    1: normal toric variety
    2: normal toric variety
    3: normal toric variety
  in the coordinate(s)
    1: [x_1_1, x_2_1]
    2: [x_1_2, x_2_2]
    3: [x_1_3, x_2_3]
```
"""
@attr Any function underlying_scheme(Z::NormalToricVariety)
  @req is_pure(polyhedral_fan(Z)) "underlying_scheme is currently only supported for toric varieties whose fan is pure"
  patch_list = affine_open_covering(Z)
  for (k, A) in enumerate(patch_list)
    C = cone(pm_object(A).WEIGHT_CONE)
    n = length(hilbert_basis(C))
    R, _ = polynomial_ring(QQ, ["x_$(i)_$(k)" for i in 1:n]; cached=false);
    set_attribute!(A, :toric_ideal, toric_ideal(R, A))
  end
  cov = Covering(patch_list)

  for i in 1:(length(patch_list)-1)
    for j in i+1:length(patch_list)
      X = patch_list[i]
      Y = patch_list[j]
      gd = ToricGluingData(Z, X, Y)
      add_gluing!(cov, LazyGluing(X, Y, gd))
      continue
      facet = intersect(cone(X), cone(Y))
      (dim(facet) == dim(cone(X)) - 1) || continue
      vmat = _find_localization_element(cone(X), cone(Y), facet)
      U = _localize_affine_toric_variety(X, vmat)
      V = _localize_affine_toric_variety(Y, (-1)*vmat)
      add_gluing!(cov, _compute_gluings(X, Y, vmat, U, V))
    end
  end

  # TODO: Improve the gluing (lazy gluing) or try to use the Hasse diagram.
  # TODO: For now, we conjecture, that the composition of the computed gluings is sufficient to deduce all gluings.
  #fill_transitions!(cov)
  return CoveredScheme(cov)
end


function _find_localization_element(X::Cone{QQFieldElem}, Y::Cone{QQFieldElem}, facet::Cone{QQFieldElem})
  CXdual = polarize(X)
  CYdual = polarize(Y)
  candidates = polarize(facet).pm_cone.HILBERT_BASIS_GENERATORS[2]
  pos = findfirst(j -> ((candidates[j,:] in CXdual) && ((-candidates[j,:]) in CYdual)), 1:nrows(candidates))
  sign = 1
  if pos === nothing
    pos = findfirst(j -> (((-candidates[j,:]) in CXdual) && (candidates[j,:] in CYdual)), 1:nrows(candidates))
    sign = -1
  end
  @req pos !== nothing "no element found for localization"
  return sign * matrix(ZZ.(collect(candidates[pos,:])))
end


function _localize_affine_toric_variety(X::AffineNormalToricVariety, vmat::ZZMatrix)
  AX = transpose(hilbert_basis(X))
  Id = identity_matrix(ZZ, ncols(AX))
  sol = solve_mixed(AX, vmat, Id, zero(vmat))
  fX = prod([x^k for (x, k) in zip(gens(ambient_coordinate_ring(X)), sol)])
  return PrincipalOpenSubset(X, OO(X)(fX))
end


function _compute_gluings(X::AffineNormalToricVariety, Y::AffineNormalToricVariety, vmat::ZZMatrix, U::PrincipalOpenSubset, V::PrincipalOpenSubset)
  AX = transpose(hilbert_basis(X))
  AY = transpose(hilbert_basis(Y))
  img_gens = _compute_image_generators(AX, AY, vmat)
  fres = hom(OO(V), OO(U), [prod([(k >= 0 ? x^k : inv(x)^(-k)) for (x, k) in zip(gens(OO(U)), w)]) for w in img_gens])
  img_gens = _compute_image_generators(AY, AX, (-1)*vmat)
  gres = hom(OO(U), OO(V), [prod([(k >= 0 ? x^k : inv(x)^(-k)) for (x, k) in zip(gens(OO(V)), w)]) for w in img_gens])
  set_attribute!(gres, :inverse, fres)
  set_attribute!(fres, :inverse, gres)
  f = morphism(U, V, fres)
  g = morphism(V, U, gres)
  set_attribute!(g, :inverse, f)
  set_attribute!(f, :inverse, g)
  return SimpleGluing(X, Y, f, g)
end


function _compute_image_generators(AX::ZZMatrix, AY::ZZMatrix, vmat::ZZMatrix)
  l = findfirst(j -> (vmat == AX[:,j]), 1:ncols(AX))
  Idext = identity_matrix(ZZ, ncols(AX))
  Idext[l,l] = 0
  img_gens = [solve_mixed(AX, AY[:, k], Idext, zero(matrix_space(ZZ, ncols(Idext), 1))) for k in 1:ncols(AY)]
end

is_irreducible(X::NormalToricVariety) = true
is_reduced(X::NormalToricVariety) = true
is_empty(X::NormalToricVariety) = false
is_integral(X::NormalToricVariety) = true
is_connected(X::NormalToricVariety) = true
