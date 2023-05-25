#######################################
### 1: Underlying scheme
#######################################

@doc raw"""
    underlying_scheme(X::ToricSpec)

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

julia> affine_toric_scheme = ToricSpec(antv)
Spec of an affine toric variety with cone spanned by RayVector{QQFieldElem}[[1, 0], [0, 1]]

julia> underlying_scheme(affine_toric_scheme)
Spec of Quotient of Multivariate polynomial ring in 2 variables over QQ by ideal(0)
```
"""
@attr underlying_scheme(X::ToricSpec) = Spec(base_ring(toric_ideal(X)), toric_ideal(X))


@doc raw"""
    underlying_scheme(X::ToricCoveredScheme)

For a toric covered scheme ``X``, this returns
the underlying scheme. In other words, by applying
this method, you obtain a scheme that has forgotten
its toric origin.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> toric_scheme = ToricCoveredScheme(P2)
Scheme of a toric variety with fan spanned by RayVector{QQFieldElem}[[1, 0], [0, 1], [-1, -1]]

julia> underlying_scheme(toric_scheme)
covered scheme with 3 affine patches in its default covering
```
"""
@attr function underlying_scheme(Z::ToricCoveredScheme)
  antv = affine_open_covering(Z)
  var_names = [Symbol.(["x_$(i)_$(k)" for i in 1:ngens(base_ring(toric_ideal(antv[k])))]) for k in 1:length(antv)]
  patch_list = [ToricSpec(U, var_names=n) for (U, n) in zip(antv, var_names)]
  cov = Covering(patch_list)
  
  for i in 1:(length(patch_list)-1)
    for j in i+1:length(patch_list)
      
      X = patch_list[i]
      Y = patch_list[j]
      
      # Step 1: Only glue if the cone of intersection has the "right" dimension:
      facet = intersect(cone(X), cone(Y))
      (dim(facet) == dim(cone(X)) - 1) || continue
      
      # Step 2: Find localization element vmat
      CXdual = dual_cone(X)
      CYdual = dual_cone(Y)
      candidates = polarize(facet).pm_cone.HILBERT_BASIS_GENERATORS[2]
      pos = findfirst(j -> ((candidates[j,:] in CXdual) && ((-candidates[j,:]) in CYdual)), 1:nrows(candidates))
      sign = 1
      if pos === nothing
        pos = findfirst(j -> (((-candidates[j,:]) in CXdual) && (candidates[j,:] in CYdual)), 1:nrows(candidates))
        sign = -1
      end
      @req pos !== nothing "no element found for localization"
      vmat = sign * matrix(ZZ.(collect(candidates[pos,:])))
      
      # Step 3: Localize X at vmat
      AX = transpose(hilbert_basis(X))
      Id = identity_matrix(ZZ, ncols(AX))
      sol = solve_mixed(AX, vmat, Id, zero(vmat))
      fX = prod([x^k for (x, k) in zip(gens(ambient_coordinate_ring(X)), sol)])
      U = PrincipalOpenSubset(X, OO(X)(fX))
      
      # Step 4: Localize Y at vmat
      AY = transpose(hilbert_basis(Y))
      Id = identity_matrix(ZZ, ncols(AY))
      sol = solve_mixed(AY, -vmat, Id, zero(vmat))
      fY = prod([x^(k>0 ? 1 : 0) for (x, k) in zip(gens(ambient_coordinate_ring(Y)), sol)])
      V = PrincipalOpenSubset(Y, OO(Y)(fY))
      
      # Step 5: Compute the glueing isomorphisms V -> U
      l = findfirst(j -> (vmat == AX[:,j]), 1:ncols(AX))
      Idext = identity_matrix(ZZ, ncols(AX))
      Idext[l,l] = 0
      img_gens = [solve_mixed(AX, AY[:, k], Idext, zero(matrix_space(ZZ, ncols(Idext), 1))) for k in 1:ncols(AY)]
      fres = hom(OO(V), OO(U), [prod([(k >= 0 ? x^k : inv(x)^(-k)) for (x, k) in zip(gens(OO(U)), w)]) for w in img_gens])
      
      # Step 6: Compute the glueing isomorphisms U -> V
      l = findfirst(j -> (-vmat == AY[:,j]), 1:ncols(AY))
      Idext = identity_matrix(ZZ, ncols(AY))
      Idext[l,l] = 0
      img_gens = [solve_mixed(AY, AX[:, k], Idext, zero(matrix_space(ZZ, ncols(Idext), 1))) for k in 1:ncols(AX)]
      gres = hom(OO(U), OO(V), [prod([(k >= 0 ? x^k : inv(x)^(-k)) for (x, k) in zip(gens(OO(V)), w)]) for w in img_gens])
      
      # Step 7: Assign the gluing isomorphisms
      set_attribute!(gres, :inverse, fres)
      set_attribute!(fres, :inverse, gres)
      f = SpecMor(U, V, fres)
      g = SpecMor(V, U, gres)
      set_attribute!(g, :inverse, f)
      set_attribute!(f, :inverse, g)
      G = SimpleGlueing(X, Y, f, g)
      add_glueing!(cov, G)
      
    end
  end
  
  # TODO: Improve the gluing (lazy gluing) or try to use the Hasse diagram.
  # TODO: For now, we conjecture, that the composition of the computed glueings is sufficient to deduce all glueings.
  fill_transitions!(cov)
  return CoveredScheme(cov)
end



#######################################
### 2: Attributes of ToricCoveredScheme
#######################################

affine_open_covering(X::ToricCoveredScheme) = affine_open_covering(X.ntv)
character_lattice(X::ToricCoveredScheme) = character_lattice(X.ntv)
character_to_rational_function(R::MPolyRing, X::ToricCoveredScheme, character::Vector{ZZRingElem}) = character_to_rational_function(R, X.ntv, character)
character_to_rational_function(X::ToricCoveredScheme, character::Vector{ZZRingElem}) = character_to_rational_function(X.ntv, character)
class_group(X::ToricCoveredScheme) = class_group(X.ntv)
coefficient_ring(X::ToricCoveredScheme) = coefficient_ring(X.ntv)
coordinate_names(X::ToricCoveredScheme) = coordinate_names(X.ntv)
coordinate_names_of_torus(X::ToricCoveredScheme) = coordinate_names_of_torus(X.ntv)
coordinate_ring_of_torus(R::MPolyRing, X::ToricCoveredScheme) = coordinate_ring_of_torus(R, X.ntv)
coordinate_ring_of_torus(X::ToricCoveredScheme) = coordinate_ring_of_torus(X.ntv)
cox_ring(R::MPolyRing, X::ToricCoveredScheme) = cox_ring(R, X.ntv)
cox_ring(X::ToricCoveredScheme) = cox_ring(X.ntv)
dim(X::ToricCoveredScheme) = dim(X.ntv)
dim_of_torusfactor(X::ToricCoveredScheme) = dim_of_torusfactor(X.ntv)
euler_characteristic(X::ToricCoveredScheme) = euler_characteristic(X.ntv)
fan(X::ToricCoveredScheme) = fan(X.ntv)
ideal_of_linear_relations(R::MPolyRing, X::ToricCoveredScheme) = ideal_of_linear_relations(R, X.ntv)
ideal_of_linear_relations(X::ToricCoveredScheme) = ideal_of_linear_relations(X.ntv)
irrelevant_ideal(R::MPolyRing, X::ToricCoveredScheme) = irrelevant_ideal(R, X.ntv)
irrelevant_ideal(X::ToricCoveredScheme) = irrelevant_ideal(X.ntv)
map_from_character_lattice_to_torusinvariant_weil_divisor_group(X::ToricCoveredScheme) = map_from_character_lattice_to_torusinvariant_weil_divisor_group(X.ntv)
map_from_torusinvariant_cartier_divisor_group_to_picard_group(X::ToricCoveredScheme) = map_from_torusinvariant_cartier_divisor_group_to_picard_group(X.ntv)
map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(X::ToricCoveredScheme) = map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(X.ntv)
map_from_torusinvariant_weil_divisor_group_to_class_group(X::ToricCoveredScheme) = map_from_torusinvariant_weil_divisor_group_to_class_group(X.ntv)
mori_cone(X::ToricCoveredScheme) = mori_cone(X.ntv)
nef_cone(X::ToricCoveredScheme) = nef_cone(X.ntv)
picard_group(X::ToricCoveredScheme) = picard_group(X.ntv)
set_coordinate_names(X::ToricCoveredScheme) = set_coordinate_names(X.ntv)
set_coordinate_names_of_torus(X::ToricCoveredScheme) = set_coordinate_names_of_torus(X.ntv)
stanley_reisner_ideal(R::MPolyRing, X::ToricCoveredScheme) = stanley_reisner_ideal(R, X.ntv)
stanley_reisner_ideal(X::ToricCoveredScheme) = stanley_reisner_ideal(X.ntv)
torusinvariant_cartier_divisor_group(X::ToricCoveredScheme) = torusinvariant_cartier_divisor_group(X.ntv)
torusinvariant_prime_divisors(X::ToricCoveredScheme) = torusinvariant_prime_divisors(X.ntv)
torusinvariant_weil_divisor_group(X::ToricCoveredScheme) = torusinvariant_weil_divisor_group(X.ntv)
underlying_toric_variety(X::ToricCoveredScheme) = X.ntv


#######################################
### 3: Attributes of ToricSpec
#######################################

affine_open_covering(X::ToricSpec) = affine_open_covering(X.antv)
character_lattice(X::ToricSpec) = character_lattice(X.antv)
character_to_rational_function(R::MPolyRing, X::ToricSpec, character::Vector{ZZRingElem}) = character_to_rational_function(R, X.antv, character)
character_to_rational_function(X::ToricSpec, character::Vector{ZZRingElem}) = character_to_rational_function(X.antv, character)
class_group(X::ToricSpec) = class_group(X.antv)
coefficient_ring(X::ToricSpec) = coefficient_ring(X.antv)
coordinate_names(X::ToricSpec) = coordinate_names(X.antv)
coordinate_names_of_torus(X::ToricSpec) = coordinate_names_of_torus(X.antv)
coordinate_ring_of_torus(R::MPolyRing, X::ToricSpec) = coordinate_ring_of_torus(R, X.antv)
coordinate_ring_of_torus(X::ToricSpec) = coordinate_ring_of_torus(X.antv)
cox_ring(R::MPolyRing, X::ToricSpec) = cox_ring(R, X.antv)
cox_ring(X::ToricSpec) = cox_ring(X.antv)
dim(X::ToricSpec) = dim(X.antv)
dim_of_torusfactor(X::ToricSpec) = dim_of_torusfactor(X.antv)
euler_characteristic(X::ToricSpec) = euler_characteristic(X.antv)
fan(X::ToricSpec) = fan(X.antv)
ideal_of_linear_relations(R::MPolyRing, X::ToricSpec) = ideal_of_linear_relations(R, X.antv)
ideal_of_linear_relations(X::ToricSpec) = ideal_of_linear_relations(X.antv)
irrelevant_ideal(R::MPolyRing, X::ToricSpec) = irrelevant_ideal(R, X.antv)
irrelevant_ideal(X::ToricSpec) = irrelevant_ideal(X.antv)
map_from_character_lattice_to_torusinvariant_weil_divisor_group(X::ToricSpec) = map_from_character_lattice_to_torusinvariant_weil_divisor_group(X.antv)
map_from_torusinvariant_cartier_divisor_group_to_picard_group(X::ToricSpec) = map_from_torusinvariant_cartier_divisor_group_to_picard_group(X.antv)
map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(X::ToricSpec) = map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(X.antv)
map_from_torusinvariant_weil_divisor_group_to_class_group(X::ToricSpec) = map_from_torusinvariant_weil_divisor_group_to_class_group(X.antv)
mori_cone(X::ToricSpec) = mori_cone(X.antv)
nef_cone(X::ToricSpec) = nef_cone(X.antv)
picard_group(X::ToricSpec) = picard_group(X.antv)
set_coordinate_names(X::ToricSpec) = set_coordinate_names(X.antv)
set_coordinate_names_of_torus(X::ToricSpec) = set_coordinate_names_of_torus(X.antv)
stanley_reisner_ideal(R::MPolyRing, X::ToricSpec) = stanley_reisner_ideal(R, X.antv)
stanley_reisner_ideal(X::ToricSpec) = stanley_reisner_ideal(X.antv)
torusinvariant_cartier_divisor_group(X::ToricSpec) = torusinvariant_cartier_divisor_group(X.antv)
torusinvariant_prime_divisors(X::ToricSpec) = torusinvariant_prime_divisors(X.antv)
torusinvariant_weil_divisor_group(X::ToricSpec) = torusinvariant_weil_divisor_group(X.antv)
underlying_toric_variety(X::ToricSpec) = X.antv

# Attributes special to affine toric varieties
cone(X::ToricSpec) = cone(X.antv)
dual_cone(X::ToricSpec) = dual_cone(X.antv)
hilbert_basis(X::ToricSpec) = hilbert_basis(X.antv)
toric_ideal(R::MPolyRing, X::ToricSpec) = toric_ideal(R, X.antv)
toric_ideal(X::ToricSpec) = toric_ideal(X.antv)
