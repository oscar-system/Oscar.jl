@testset "determinants of complexes" begin
  R, (x, y, z) = QQ[:x, :y, :z]
  I = ideal(R, gens(R))
  k, _ = quo(R, I)
  M = quotient_ring_as_module(k)
  res, _ = free_resolution(Oscar.SimpleFreeResolution, M)
  p = det(res)
  @test is_one(p)
  @test p == det(res; direction=:from_right_to_left, upper_bound=5)
  c = hom(res, free_module(R, 1))
  @test p == det(c; lower_bound=-5);
  @test p == det(c; direction=:from_right_to_left)
end

