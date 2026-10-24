# more tests are in `test/AlgebraicGeometry/Schemes/`
# and in `/test/AlgebraicGeometry/ToricVarieties/`
#
@testset "Blowups in existing rays" begin
  c = positive_hull([1 0 0; 1 1 0; 1 1 1; 1 0 1])
  ntv = normal_toric_variety(c)
  bl0 = blow_up(ntv, [1,0,0])
  bl1 = blow_up(domain(bl0), [1,0,0])
  @test n_maximal_cones(ntv) == 1
  @test n_maximal_cones(domain(bl0)) == 2
  @test n_maximal_cones(domain(bl1)) == 2
  @test n_rays(ntv) == 4
  @test n_rays(domain(bl0)) == 4
  @test n_rays(domain(bl1)) == 4
  @test exceptional_prime_divisor(bl0) == toric_divisor(domain(bl0), [1,0,0,0])
  @test exceptional_prime_divisor(bl1) == toric_divisor(domain(bl1), [1,0,0,0])
end

@testset "Blowups along minimal supercone coordinates" begin
  P2 = projective_space(NormalToricVariety, 2)
  f_1 = blow_up_along_minimal_supercone_coordinates(P2, [2, 3, 0])
  f_2 = blow_up(P2, [2, 3])
  pf_1 = polyhedral_fan(domain(f_1))
  pf_2 = polyhedral_fan(domain(f_2))
  @test Set(rays(pf_1)) == Set(rays(pf_2))
end

@testset "Total and strict transforms in Cox rings" begin
  # Affine space, (2, 3)-blowup
  X = affine_space(NormalToricVariety, 2)
  set_coordinate_names(X, [:x, :y])
  R = cox_ring(X)
  x, y = gens(R)
  f = blow_up(X, [2, 3]; coordinate_name=:u)
  Y = domain(f)
  S = cox_ring(Y)
  x_, y_, u = gens(S)

  ## Polynomial test
  g = x + y^3
  @test cox_ring_module_homomorphism(f, g) == x_*u^2 + y_^3*u^9

  ## Subscheme is a curve
  I = ideal(R, [x^3 + y^3])
  J, k = strict_transform_with_index(f, I)
  @test J == ideal(S, [x_^3 + y_^3*u^3])
  @test k == 6
  g = x^3 + y^3
  h, k_poly = strict_transform_with_index(f, g)
  @test ideal(h) == J
  @test k == k_poly

  # 1/2(1, 1) quotient singularity, (1/2, 1/2)-blowup
  ray_generators = [[2, -1], [0, 1]]
  max_cones = IncidenceMatrix([[1, 2]])
  X = normal_toric_variety(max_cones, ray_generators)
  set_coordinate_names(X, ["x", "y"])
  R = cox_ring(X)
  x, y = gens(R)
  f = blow_up(X, [1, 0]; coordinate_name=:u)
  Y = domain(f)
  S = cox_ring(Y)
  x_, y_, u = gens(S)

  ## Polynomial test
  g = x + y^3
  @test cox_ring_module_homomorphism(f, g) == x_*u + y_^3*u^2

  ## Subscheme is a curve
  I = ideal(R, [x + y^3])
  J_strict = strict_transform(f, I)
  @test J_strict == ideal(S, [x_ + y_^3*u])
  @test ideal_sheaf(Y, J_strict) == strict_transform(f, ideal_sheaf(X, I))
  J_total = total_transform(f, I)
  @test J_total == ideal(S, [x_*u + y_^3*u^2])
  @test ideal_sheaf(Y, J_total) == total_transform(f, ideal_sheaf(X, I))
  g = x + y^3
  h_strict = strict_transform(f, g)
  @test ideal(h_strict) == J_strict
  h_total = total_transform(f, g)
  @test ideal(h_total) == J_total

  ## Subscheme is a zero-dimensional scheme, topologically three points
  I = ideal(R, [x - y^3, x - y^5])
  J_strict = strict_transform(f, I)
  @test J_strict == ideal(S, [x_^2 - x_*y_, x_*y_^2 - y_^3, -x_ + y_^3*u])
  @test ideal_sheaf(Y, J_strict) == strict_transform(f, ideal_sheaf(X, I))
  J_total = total_transform(f, I)
  @test J_total == ideal(S, [x_*u - y_^3*u^2, x_*u - y_^5*u^3])
  @test ideal_sheaf(Y, J_total) == total_transform(f, ideal_sheaf(X, I))

  # Quadratic cone, blowup along origin
  ray_generators = [[0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]
  max_cones = IncidenceMatrix([[1, 2, 3, 4]])
  PF = polyhedral_fan(max_cones, ray_generators)
  X = normal_toric_variety(PF)
  set_coordinate_names(X, [:x, :y, :z, :w])
  R = cox_ring(X)
  x, y, z, w = gens(R)
  f = blow_up(X, [1, 1, 2]; coordinate_name=:u)
  Y = domain(f)
  S = cox_ring(Y)
  x_, y_, z_, w_, u = gens(S)

  ## Subscheme is a curve
  I = ideal(R, [x + z^2*y])
  J_strict = strict_transform(f, I)
  @test J_strict == ideal(S, [x_ + y_*z_^2*u])
  @test ideal_sheaf(Y, J_strict) == strict_transform(f, ideal_sheaf(X, I))

  # Quadratic cone, blowup along non-Cartier divisor
  ray_generators = [[0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]
  max_cones = IncidenceMatrix([[1, 2, 3, 4]])
  PF = polyhedral_fan(max_cones, ray_generators)
  X = normal_toric_variety(PF)
  set_coordinate_names(X, [:x, :y, :z, :w])
  R = cox_ring(X)
  x, y, z, w = gens(R)
  f = blow_up_along_minimal_supercone_coordinates(X, [1, 0, 0, 0]; coordinate_name=:u)
  Y = domain(f)
  S = cox_ring(Y)
  x_, y_, z_, w_ = gens(S)

  ## Subscheme is a curve
  I = ideal(R, [x + z^2*y])
  J_strict = strict_transform(f, I)
  @test J_strict == ideal(S, [x_ + y_*z_^2])
  @test ideal_sheaf(Y, J_strict) == strict_transform(f, ideal_sheaf(X, I))
end

