@testset "basics about blowups" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  f = x^2 + y^3 + z^5
  X = spec(R, ideal(R, f))
  I = ideal(OO(X), [x, y, z])
  IZ = IdealSheaf(X, I)

  bl = blow_up(IZ)

  @test bl isa AbsCoveredSchemeMorphism{<:AbsCoveredScheme, typeof(space(IZ)), Nothing, <:BlowupMorphism}
  @test Oscar.underlying_morphism(bl) === projection(bl)

  Y = domain(bl)
  @test codomain(bl) === space(IZ)
  @test Y isa AbsCoveredScheme

  E = exceptional_divisor(bl)
end

@testset "strict transforms of cartier divisors" begin
  IP2 = projective_space(QQ, [:x, :y, :z])
  S = ambient_coordinate_ring(IP2)
  (x,y,z) = gens(S)
  I = ideal(S, [x, y])
  #set_name!(IP2, "ℙ²")
  II = IdealSheaf(IP2, I)
  p = blow_up(II)
  C = Oscar.effective_cartier_divisor(IP2, (x+y)^2)
  D = Oscar.effective_cartier_divisor(IP2, (x-3*y)^5)
  C_trans = strict_transform(p, C)
  @test weil_divisor(pullback(projection(p))(C)) - weil_divisor(C_trans) == weil_divisor(2*exceptional_divisor(p))
  D_trans = strict_transform(p, D)
  @test weil_divisor(pullback(projection(p))(D)) - weil_divisor(D_trans) == weil_divisor(5*exceptional_divisor(p))
  CD = 5*C + 7*D
  CD_trans = strict_transform(p, CD)
  @test weil_divisor(CD_trans) == weil_divisor(5*C_trans + 7*D_trans)
end

@testset "isomorphism on complement of center" begin
  P = projective_space(QQ, [:x, :y, :z])
  S = homogeneous_coordinate_ring(P)
  (x, y, z) = gens(S)
  II = IdealSheaf(P, [x, y])
  JJ = IdealSheaf(P, [x^2*z-y^3])
  Y = scheme(II)
  p = blow_up(II)
  X = domain(p)
  f = Oscar.isomorphism_on_complement_of_center(p)
  h = inverse(f)
  U = domain(f)
  @test compose(f, h) == identity_map(U)
  V = codomain(f)
  @test compose(h, f) == identity_map(V)
  KY = function_field(Y)
  KX = function_field(X)
  y, z = gens(ambient_coordinate_ring(first(affine_charts(Y))))
  a = KY(y, z)
  b = KY(a[affine_charts(Y)[2]])
  @test pullback(p)(a) == pullback(p)(b)

  inc_C = Oscar.CoveredClosedEmbedding(scheme(JJ), JJ)
  C = domain(inc_C)
  C_up, inc_C_up, p_res = strict_transform(p, inc_C)


  KC_up = function_field(C_up)
  KC = function_field(C)
  aa = KC(a[affine_charts(Y)[3]])
  @test pullback(p_res)(aa)^2 + one(KC_up) == pullback(p_res)(aa^2 + one(KC))

end

@testset "pushforward of function field elements through resolutions" begin
  IA3 = affine_space(QQ, [:x, :y, :z])
  x, y, z = gens(OO(IA3))
  pr1 = blow_up(IA3, ideal(OO(IA3), [x, y, z]))
  Y0 = codomain(pr1)
  Y1 = domain(pr1)
  KY0 = function_field(Y0)
  xx = KY0(x)
  yy = KY0(y)
  zz = KY0(z)

  f = x^4 - y*z
  I = ideal(OO(IA3), f)
  II = ideal_sheaf(Y0, IA3, [f])
  inc = Oscar.CoveredClosedEmbedding(Y0, II)
  X0 = domain(inc)
  KX0 = function_field(X0)
  xx0 = KX0(x)
  yy0 = KX0(y)
  zz0 = KX0(z)
  X1, inc1, pr1_res = strict_transform(pr1, inc)
  xx1 = pullback(pr1_res)(xx0)
  yy1 = pullback(pr1_res)(yy0)
  zz1 = pullback(pr1_res)(zz0)

  I_sing = radical(pushforward(inc1, Oscar.ideal_sheaf_of_singular_locus(X1)))

  pr2 = blow_up(I_sing)
  @test scheme(I_sing) === domain(pr1)
  @test codomain(pr2) === domain(pr1)
  Y2 = domain(pr2)
  X2, inc2, pr2_res = strict_transform(pr2, inc1)

  pr_res = Oscar.composite_map(pr2_res, pr1_res)
  pr = Oscar.compose(pr2, pr1)

  pr_inc = compose(compose(inc2, pr2), pr1)
  @test Oscar.underlying_morphism(pr1) isa CoveredSchemeMorphism
  @test domain(Oscar.underlying_morphism(pr1)) === domain(pr1)
  @test codomain(Oscar.underlying_morphism(pr1)) === codomain(pr1)
  @test compose(inc2, pr2) isa Oscar.CompositeCoveredSchemeMorphism
  @test domain(compose(inc2, pr2)) === domain(inc2)
  @test codomain(compose(inc2, pr2)) === codomain(pr2)

  @test domain(pr1) === codomain(compose(inc2, pr2))
  @test domain(Oscar.underlying_morphism(pr1)) === codomain(compose(inc2, pr2))
  pr_inc2 = compose(compose(inc2, pr2), Oscar.underlying_morphism(pr1))
  pr_inc3 = compose(Oscar.underlying_morphism(compose(inc2, pr2)), Oscar.underlying_morphism(pr1))
  pr_inc4 = compose(Oscar.underlying_morphism(compose(inc2, pr2)), pr1)

  @test CoveredSchemeMorphism(pr_inc) == CoveredSchemeMorphism(pr_inc2)
  @test CoveredSchemeMorphism(pr_inc) == CoveredSchemeMorphism(pr_inc3)
  @test CoveredSchemeMorphism(pr_inc) == CoveredSchemeMorphism(pr_inc4)

  pr_inc_alt = Oscar.composite_map(pr2_res, compose(inc1, pr1))
  pr_inc_alt2 = Oscar.composite_map(compose(pr2_res, inc1), pr1)
  pr_inc_alt3 = Oscar.composite_map(pr2_res, compose(inc1, Oscar.underlying_morphism(pr1)))
  pr_inc_alt4 = Oscar.composite_map(compose(pr2_res, inc1), Oscar.underlying_morphism(pr1))

  @test CoveredSchemeMorphism(pr_inc) == CoveredSchemeMorphism(pr_inc_alt)
  @test CoveredSchemeMorphism(pr_inc) == CoveredSchemeMorphism(pr_inc_alt2)
  @test CoveredSchemeMorphism(pr_inc) == CoveredSchemeMorphism(pr_inc_alt3)
  @test CoveredSchemeMorphism(pr_inc) == CoveredSchemeMorphism(pr_inc_alt4)
end

