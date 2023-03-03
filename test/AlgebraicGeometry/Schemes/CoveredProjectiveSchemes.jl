@testset "blowups" begin
  P = projective_space(QQ, ["x", "y", "z", "w"])
  X = covered_scheme(P)
  S = ambient_coordinate_ring(P)
  (x, y, z, w) = gens(S)
  M = S[x y z; y z w]
  I = ideal(S, minors(M, 2))
  II = IdealSheaf(P, I)
  B = blow_up(II)::BlowupMorphism
  Bcov = domain(B)
  p = projection(B)
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
  @test B1 isa AbsProjectiveScheme
  Y = covered_scheme(B1)
  @test !is_smooth(Y)
  S, inc = singular_locus(Y);
  @test dim(S) == 1
  U = affine_charts(S)
  @test !is_empty(U[1])
  @test dim(U[1]) == 1
  @test dim(U[2]) == 0
  B2 = blow_up(image_ideal(inc), var_name="t") # Use a different letter for the homogeneous variables in the 2nd blowup
  Z = domain(B2)
  @test is_smooth(Z)
end

@testset "K3 surface reconstructed" begin 
  IP1 = projective_space(GF(29), ["s", "t"])

  O1 = twisting_sheaf(IP1, -1)
  O4 = twisting_sheaf(IP1, -4)
  O6 = twisting_sheaf(IP1, -6)

  E = direct_sum([O4, O6, O1])

  X_proj = projectivization(E, var_names=["x", "y", "z"])

  X = covered_scheme(X_proj)

  U = affine_charts(X)[3]
  (x, y, t) = gens(OO(U))
  ft = y^2 - (x^3 + 21*x + (28*t^7+18))
  I = IdealSheaf(X, U, [ft])

  # The fabulous K3-surface. Almost.
  S = subscheme(I)

  sing_S, inc_sing = singular_locus(S)

  I_sing = image_ideal(inc_sing)

  @test scheme(I_sing) === S

  Y1_proj = blow_up(I_sing)
  Y1 = domain(Y1_proj)
  simplify!(Y1)
  sing_Y1, inc_Y1 = singular_locus(Y1)
  I_sing_Y1 = image_ideal(inc_Y1)
  @test !is_smooth(Y1)
end
