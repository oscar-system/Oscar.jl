# more tests are in `test/AlgebraicGeometry/Schemes/`
#
@testset "Blow ups in existing rays" begin
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

@testset "Total and strict transforms in Cox rings" begin
  # 1/2(1, 1) quotient singularity, (1/2, 1/2)-blowup
  ray_generators = [[2, -1], [0, 1]]
  max_cones = IncidenceMatrix([[1, 2]])
  X = normal_toric_variety(max_cones, ray_generators)
  set_coordinate_names(X, ["x", "y"])
  R = cox_ring(X)
  x, y = gens(R)
  f = blow_up(X, [1, 0]; coordinate_name="u")
  Y = domain(f)
  S = cox_ring(Y)
  x_, y_, u = gens(S)
  I = ideal(R, [x + y^3])
  @test strict_transform_with_index(f, I) == (
    ideal(S, [x_ + y_^3*u]), 1
  )
  @test (
    ideal_sheaf(domain(f), strict_transform(f, I))
    == strict_transform(f, ideal_sheaf(codomain(f), I))
  )
  I = ideal(R, [x + y^3, x - y^5])
  @test strict_transform_with_index(f, I) == (
    ideal(S, [x_^2 - x_*y_, x_*y_^2 - y_^3, x_ + y_^3*u]), 2
  )
  @test (
    ideal_sheaf(domain(f), strict_transform(f, I))
    == strict_transform(f, ideal_sheaf(codomain(f), I))
  )
end
