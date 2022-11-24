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

@testset "blowups II" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  A3 = Spec(R)
  M = R[x y z; y-1 z-2 x-3]
  I = ideal(R, minors(M, 2))
  BlA3 = blow_up(A3, I)
  p = covered_projection_to_base(BlA3)
  @test domain(p) === covered_scheme(BlA3)
  @test codomain(p)[1][1] === A3
  @test is_smooth(covered_scheme(BlA3))

  J = ideal(R, [x,y,z])
  Blo = blow_up(A3, J)
  p2 = covered_projection_to_base(Blo)
  @test domain(p2) === covered_scheme(Blo)
  @test codomain(p2)[1][1] === A3
  @test is_smooth(covered_scheme(Blo))
end

