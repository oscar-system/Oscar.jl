@testset "degree zero complexes" begin
  S, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
  F = graded_free_module(S, [0, 0, 0])
  v = x*F[1] + y*F[2] + z*F[3]
  kosz = koszul_complex(Oscar.KoszulComplex, v)
  c = Oscar.DegreeZeroComplex(kosz)
  @test c[0] isa FreeMod
  @test !Oscar.can_compute_index(c, -1)
  @test Oscar.can_compute_index(c, 3)
  phi = map(c, 2)
  @test is_homogeneous(phi)
  phi = map(c, 3)
  @test is_homogeneous(phi)
  @test Oscar.has_upper_bound(c)
end
