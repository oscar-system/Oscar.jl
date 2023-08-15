@testset "basics about blowups" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  f = x^2 + y^3 + z^5
  X = CoveredScheme(Spec(R, ideal(R, f)))
  U = X[1][1] # the first chart

  IZ = IdealSheaf(X, U, OO(U).([x, y, z]))

  bl = blow_up(IZ)

  @test bl isa AbsCoveredSchemeMorphism{<:AbsCoveredScheme, typeof(X), Nothing, BlowupMorphism}
  @test underlying_morphism(bl) === projection(bl)

  Y = domain(bl)
  @test codomain(bl) === X
  @test Y isa AbsCoveredScheme

  E = exceptional_divisor(bl)
end

@testset "strict transforms of cartier divisors" begin
  IP2 = projective_space(QQ, ["x", "y", "z"])
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
  P = projective_space(QQ, ["x", "y", "z"])
  S = homogeneous_coordinate_ring(P)
  (x, y, z) = gens(S)
  II = IdealSheaf(P, [x, y])
  Y = scheme(II)
  p = blow_up(II)
  X = domain(p)
  f = oscar.isomorphism_on_complement_of_center(p)
  h = inverse(f)
  U = domain(f)
  @test compose(f, h) == identity_map(U)
  V = codomain(f)
  @test compose(h, f) == identity_map(V)
  g = oscar.isomorphism_on_open_subset(p)
  @test is_isomorphism(g)
  KY = function_field(Y)
  KX = function_field(X)
  y, z = gens(ambient_coordinate_ring(first(affine_charts(Y))))
  a = KY(y, z)
  b = KY(a[affine_charts(Y)[2]])
  @test pullback(p)(a) == pullback(p)(b)
end

