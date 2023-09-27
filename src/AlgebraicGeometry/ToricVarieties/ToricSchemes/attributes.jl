#######################################
### Forget toric structure
#######################################

@doc raw"""
    forget_toric_structure(X::AffineNormalToricVariety)

Returns a pair `(Y, iso)` where `Y` is a scheme without toric structure, 
together with an isomorphism `iso : Y → X`.

# Examples
```jldoctest
julia> C = positive_hull([1 0; 0 1])
Polyhedral cone in ambient dimension 2

julia> antv = affine_normal_toric_variety(C)
Normal, affine toric variety

julia> forget_toric_structure(antv)
(V(0), Morphism: V(0) -> Normal, affine toric variety)
```
"""
function forget_toric_structure(X::AffineNormalToricVariety)
  Y = underlying_scheme(X)
  iso = SpecMor(Y, X, identity_map(OO(X)), check=true)
  iso_inv = SpecMor(X, Y, identity_map(OO(X)), check=true)
  set_attribute!(iso, :inverse => iso_inv)
  set_attribute!(iso_inv, :inverse => iso)
  return Y, iso
end

@doc raw"""
    forget_toric_structure(X::NormalToricVariety)

Returns a pair `(Y, iso)` where `Y` is a scheme without toric structure, 
together with an isomorphism `iso : Y → X`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> forget_toric_structure(P2)
(Scheme over QQ covered with 3 patches, Morphism: scheme over QQ covered with 3 patches -> normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor)
```
"""
function forget_toric_structure(X::NormalToricVariety)
  # Collect all the isomorphisms forgetting the toric structure
  iso_dict = IdDict{AbsSpec, AbsSpecMor}()
  for U in affine_charts(X)
    iso_dict[U] = forget_toric_structure(U)[2] # store only the isomorphism
  end
  cov = Covering([domain(phi) for (U, phi) in iso_dict])

  # Prepare a dictionary that can be used in the constructor of the covering morphism
  iso_dict_covariant = IdDict{AbsSpec, AbsSpecMor}()
  for (U, phi) in iso_dict
    iso_dict_covariant[domain(phi)] = phi
  end

  # Recreate all the glueings with the new patches
  for U in affine_charts(X), V in affine_charts(X)
    glue = default_covering(X)[U, V]
    new_glue = restrict(glue, inverse(iso_dict[U]), inverse(iso_dict[V]), check=true)
    add_glueing!(cov, new_glue)
  end

  # Prepare the underlying covering morphisms for the identifying isomorphisms
  iso_cov = CoveringMorphism(cov, default_covering(X), iso_dict_covariant)
  inv_dict = IdDict{AbsSpec, AbsSpecMor}()
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
Normal, affine toric variety

julia> underlying_scheme(antv)
Spectrum
  of quotient
    of multivariate polynomial ring in 2 variables x1, x2
      over rational field
    by ideal(0)
```
"""
@attr Spec{QQField, MPolyQuoRing{QQMPolyRingElem}} underlying_scheme(X::AffineNormalToricVariety) = Spec(base_ring(toric_ideal(X)), toric_ideal(X))

### 
# Some additional structure to make computation of toric glueings lazy
struct ToricGlueingData
  X::NormalToricVariety
  U::AffineNormalToricVariety
  V::AffineNormalToricVariety
# i::Int # The indices of the charts to be glued
# j::Int
end

function _compute_toric_glueing(gd::ToricGlueingData)
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
  # elements of these Hilbert bases. The glueing isomorphisms 
  # are then obtained by expressing the generators on the one 
  # side in terms of the others. 

  hb_L = tau_dual.pm_cone.HILBERT_BASIS_GENERATORS[2] 
  
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
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> underlying_scheme(P2)
Scheme
  over rational field
with default covering
  described by patches
    1: normal, affine toric variety
    2: normal, affine toric variety
    3: normal, affine toric variety
  in the coordinate(s)
    1: [x_1_1, x_2_1]
    2: [x_1_2, x_2_2]
    3: [x_1_3, x_2_3]
```
"""
@attr function underlying_scheme(Z::NormalToricVariety)
  @req is_pure(polyhedral_fan(Z)) "underlying_scheme is currently only supported for toric varieties whose fan is pure"
  patch_list = affine_open_covering(Z)
  for (k, A) in enumerate(patch_list)
    C = cone(pm_object(A).WEIGHT_CONE)
    n = length(hilbert_basis(C))
    R, _ = polynomial_ring(QQ, ["x_$(i)_$(k)" for i in 1:n], cached = false);
    set_attribute!(A, :toric_ideal, toric_ideal(R, A))
  end
  cov = Covering(patch_list)
  
  for i in 1:(length(patch_list)-1)
    for j in i+1:length(patch_list)
      X = patch_list[i]
      Y = patch_list[j]
      gd = ToricGlueingData(Z, X, Y)
      add_glueing!(cov, LazyGlueing(X, Y, _compute_toric_glueing, gd))
      continue
      facet = intersect(cone(X), cone(Y))
      (dim(facet) == dim(cone(X)) - 1) || continue
      vmat = _find_localization_element(cone(X), cone(Y), facet)
      U = _localize_affine_toric_variety(X, vmat)
      V = _localize_affine_toric_variety(Y, (-1)*vmat)
      add_glueing!(cov, _compute_gluings(X, Y, vmat, U, V))
    end
  end
  
  # TODO: Improve the gluing (lazy gluing) or try to use the Hasse diagram.
  # TODO: For now, we conjecture, that the composition of the computed glueings is sufficient to deduce all glueings.
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
  f = SpecMor(U, V, fres)
  g = SpecMor(V, U, gres)
  set_attribute!(g, :inverse, f)
  set_attribute!(f, :inverse, g)
  return SimpleGlueing(X, Y, f, g)
end


function _compute_image_generators(AX::ZZMatrix, AY::ZZMatrix, vmat::ZZMatrix)
  l = findfirst(j -> (vmat == AX[:,j]), 1:ncols(AX))
  Idext = identity_matrix(ZZ, ncols(AX))
  Idext[l,l] = 0
  img_gens = [solve_mixed(AX, AY[:, k], Idext, zero(matrix_space(ZZ, ncols(Idext), 1))) for k in 1:ncols(AY)]
end
