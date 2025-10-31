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
  
  I = ideal(R, [x])
  k, _ = quo(R, I)
  M = quotient_ring_as_module(k)
  res, _ = free_resolution(Oscar.SimpleFreeResolution, M)
  p = det(res)
  @test p == det(res; direction=:from_right_to_left, upper_bound=5)
  c = hom(res, free_module(R, 1))
  @test p == 1//det(c; lower_bound=-5);
  @test p == 1//det(c; direction=:from_right_to_left)
  cc = Oscar.ReflectedComplex(c)
  @test p == 1//det(cc; direction=:from_right_to_left, upper_bound=5)
  @test p == 1//det(cc; direction=:from_left_to_right, lower_bound=-5)
  ccc = hom(cc, free_module(R, 1))
  @test p == det(ccc; lower_bound=-5);
  @test p == det(ccc; direction=:from_right_to_left)
end

