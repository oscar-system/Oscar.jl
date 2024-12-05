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

