@testset "blowups" begin
  P = projective_space(QQ, ["x", "y", "z", "w"])
  X = covered_scheme(P)
  S = ambient_coordinate_ring(P)
  (x, y, z, w) = gens(S)
  M = S[x y z; y z w]
  I = ideal(S, minors(M, 2))
  II = IdealSheaf(P, I)
  B = blow_up(II)
  Bcov = covered_scheme(B)
  p = covered_projection_to_base(B)
  @test domain(p) === Bcov
  @test codomain(p) === covered_scheme(P)
end
