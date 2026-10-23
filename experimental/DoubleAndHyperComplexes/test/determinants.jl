@testset "determinants of complexes" begin
  R, (x, y, z) = QQ[:x, :y, :z]
  I = ideal(R, gens(R))
  k, _ = quo(R, I)
  M = quotient_ring_as_module(k)
  res, _ = free_resolution(Oscar.SimpleFreeResolution, M)
  p = det(res; direction=:from_left_to_right)
  @test is_one(denominator(p)) && is_unit(numerator(p))
  pp = det(res; direction=:from_right_to_left, upper_bound=5)
  @test is_one(denominator(p)) && is_associated(numerator(p), numerator(pp))
  ref = Oscar.ReflectedComplex(res)
  p = det(ref; direction=:from_left_to_right, lower_bound=-5)
  pp = det(ref; direction=:from_right_to_left)
  @test is_one(denominator(p)) && is_unit(numerator(p))
  @test is_one(denominator(p)) && is_associated(numerator(p), numerator(pp))
  
  p = det(res; direction=:from_right_to_left, upper_bound=5)
  @test is_one(denominator(p)) && is_unit(numerator(p))
  c = hom(res, free_module(R, 1))
  p = det(c; lower_bound=-5);
  @test is_one(denominator(p)) && is_unit(numerator(p))
  p = det(c; direction=:from_right_to_left)
  @test is_one(denominator(p)) && is_unit(numerator(p))
  
  I = ideal(R, [x])
  k, _ = quo(R, I)
  M = quotient_ring_as_module(k)
  res, _ = free_resolution(Oscar.SimpleFreeResolution, M)
  p = det(res)
  pp = det(res; direction=:from_right_to_left, upper_bound=5)
  @test is_associated(numerator(p), numerator(pp))
  c = hom(res, free_module(R, 1))
  q = p*det(c; lower_bound=-5);
  @test is_one(denominator(q)) && is_unit(numerator(q))
  q = p*det(c; direction=:from_right_to_left)
  @test is_one(denominator(q)) && is_unit(numerator(q))
  cc = Oscar.ReflectedComplex(c)
  q = p*det(cc; direction=:from_right_to_left, upper_bound=5)
  @test is_one(denominator(q)) && is_unit(numerator(q))
  q = p*det(cc; direction=:from_left_to_right, lower_bound=0)
  @test is_one(denominator(q)) && is_unit(numerator(q))
  ccc = hom(cc, free_module(R, 1))
  q = p//det(ccc; lower_bound=-5);
  @test is_one(denominator(q)) && is_unit(numerator(q))
  q = p//det(ccc; direction=:from_right_to_left)
  @test is_one(denominator(q)) && is_unit(numerator(q))
end

