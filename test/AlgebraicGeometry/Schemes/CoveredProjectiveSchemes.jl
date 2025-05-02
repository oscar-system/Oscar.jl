@testset "blowups" begin
  P = projective_space(QQ, ["x", "y", "z", "w"])
  X = covered_scheme(P)
  S = homogeneous_coordinate_ring(P)
  (x, y, z, w) = gens(S)
  M = S[x y z; y z w]
  I = ideal(S, minors(M, 2))
  II = IdealSheaf(P, I)
  B = blow_up(II)::BlowupMorphism
  Bcov = domain(B)
  p = projection(B)
  @test domain(p) === Bcov
  @test codomain(p) === covered_scheme(P)

  X = affine_space(QQ,2)
  x1, x2 = coordinates(X)
  I = ideal([x1^2, x2^3])
  p = blow_up(X, I)
  @test !is_normal(domain(p))
end

@testset "Oscar.blow_up_chart" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  A3 = spec(R)
  M = R[x y z; y-1 z-2 x-3]
  I = ideal(R, minors(M, 2))
  BlA3 = Oscar.blow_up_chart(A3, I)
  p = covered_projection_to_base(BlA3)
  @test domain(p) === covered_scheme(BlA3)
  @test codomain(p)[1][1] === A3
  @test is_smooth(covered_scheme(BlA3))

  J = ideal(R, [x,y,z])
  Blo = Oscar.blow_up_chart(A3, J)
  p2 = covered_projection_to_base(Blo)
  @test domain(p2) === covered_scheme(Blo)
  @test codomain(p2)[1][1] === A3
  @test is_smooth(covered_scheme(Blo))
end

@testset "winter school presentation" begin
  P, (x,y,z) = QQ[:x, :y, :z]
  IA3 = spec(P)
  f = x^2-y*z^2
  I = ideal(P, f)
  X = subscheme(IA3, I)
  S, inc = singular_locus(X);
  @test S isa AbsAffineScheme
  @test inc isa ClosedEmbedding
  B1 = Oscar.blow_up_chart(X, ideal(OO(X), [x,y,z]))
  @test B1 isa AbsProjectiveScheme
  Y = covered_scheme(B1)
  @test !is_smooth(Y)
  S, inc = singular_locus(Y);
  @test dim(S) == 1
  U = affine_charts(S)
  @test !is_empty(U[1])
  @test dim(U[1]) == 1
  @test dim(U[2]) == 0
  B2 = blow_up(simplify(image_ideal(inc)), var_name="t") # Use a different letter for the homogeneous variables in the 2nd blowup
  Z = domain(B2)
  @test is_smooth(Z)
end

@testset "K3 surface reconstructed" begin
  IP1 = projective_space(GF(29), ["s", "t"])

  O1 = twisting_sheaf(IP1, -1)
  O4 = twisting_sheaf(IP1, -4)
  O6 = twisting_sheaf(IP1, -6)

  E = direct_sum([O4, O6, O1])

  X_proj = projectivization(E, var_names=[:x, :y, :z])

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

@testset "decomposition information" begin
  P3 = projective_space(QQ, 3)
  X = covered_scheme(P3)
  S = homogeneous_coordinate_ring(P3)
  (x,y, z, w) = gens(S)
  A = S[x y z; y z w]
  I = ideal(S, minors(A, 2))
  II = ideal_sheaf(P3, I)
  pr = blow_up(II)
  Y = domain(projection(pr))
  @test Oscar.has_decomposition_info(Oscar.simplified_covering(Y))
  E = exceptional_divisor(pr)
  IE = ideal_sheaf(E)
  Z = subscheme(IE)
  @test Oscar.has_decomposition_info(Oscar.simplified_covering(Z))
end
