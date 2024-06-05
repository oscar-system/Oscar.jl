@testset "free resolutions of ideals" begin
  X = veronese();
  I = defining_ideal(X);
  Pn = base_ring(I)
  FI, aug = free_resolution(Oscar.SimpleFreeResolution, I);
  @test domain(aug) === FI
  @test codomain(aug)[(0,)] isa SubquoModule
  F = graded_free_module(Pn, 1)
  dualFIC = hom(FI, F);
  # Duals invert the range from 0:10 to 0:-1:-10 etc. 
  @test dualFIC[0] isa FreeMod
  @test dualFIC[-1] isa FreeMod
  @test dualFIC[-10] isa FreeMod
  @test is_graded(dualFIC[0])
  @test is_graded(dualFIC[-1])
  @test is_graded(dualFIC[-10])
  @test !Oscar.can_compute_index(dualFIC, 1)
  betti_table(FI)
  minimal_betti_table(FI)
  simplify(FI)
end

@testset "Betti tables over quotient rings" begin
  R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
  I = ideal(R, z)
  A, _ = quo(R, I)
  (x, y, z) = gens(A)

  F = graded_free_module(A, 1)

  J, inc = sub(F, [x*F[1], y*F[1]])

  M = cokernel(inc)

  res, aug = free_resolution(Oscar.SimpleFreeResolution, M)

  betti_table(res; upper_bound=3)
  minimal_betti_table(res; upper_bound=3)
end
