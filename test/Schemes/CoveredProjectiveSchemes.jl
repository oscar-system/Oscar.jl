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

@testset "winter school presentation" begin
  P, (x,y,z) = QQ["x", "y", "z"]
  IA3 = Spec(P)
  f = x^2-y*z^2
  I = ideal(P, f)
  X = subscheme(IA3, I)
  S, inc = singular_locus(X);
  @test S isa AbsSpec
  @test inc isa ClosedEmbedding
  B1 = blow_up(X, ideal(OO(X), [x,y,z]))
  @test B1 isa ProjectiveScheme
  Y = covered_scheme(B1)
  @test !is_smooth(Y)
  S, inc = singular_locus(Y);
  @test dim(S) == 1
  U = affine_charts(S)
  @test !is_empty(U[1])
  @test dim(U[1]) == 1
  @test dim(U[2]) == 0
  B2 = blow_up(image_ideal(inc), var_name="t") # Use a different letter for the homogeneous variables in the 2nd blowup
  Z = covered_scheme(B2)
  @test is_smooth(Z)
end

