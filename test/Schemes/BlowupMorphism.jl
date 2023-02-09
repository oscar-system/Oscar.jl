@testset "basics about blowups" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  f = x^2 + y^3 + z^5
  X = CoveredScheme(Spec(R, ideal(R, f)))
  U = X[1][1] # the first chart

  IZ = IdealSheaf(X, U, OO(U).([x, y, z]))

  bl = blow_up(IZ)

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
  set_name!(X, "ℙ²")
  II = IdealSheaf(IP2, I)
  p = blow_up(II)
  C = oscar.effective_cartier_divisor(IP2, (x+y)^2)
  D = oscar.effective_cartier_divisor(IP2, (x-3*y)^5)
  C_trans, E_trans = strict_transform(p, C)
  @test E_trans == 2*exceptional_divisor(p)
  D_trans, E_trans = strict_transform(p, D)
  @test E_trans == 5*exceptional_divisor(p)
  CD = 5*C + 7*D
  CD_trans = strict_transform(p, CD)
end
