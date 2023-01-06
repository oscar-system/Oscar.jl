@testset "coherent sheaves" begin
  IP2 = projective_space(QQ, 2)
  X = covered_scheme(IP2)
  OX = StructureSheafOfRings(X)
  L = tautological_bundle(IP2)
  U = affine_charts(X)
  L(U[1])
  C = default_covering(X)
  (U12, U21) = glueing_domains(C[U[1], U[2]])
  rho = L(U[1], U21)
  @test rho(L(U[1])[1]) == inv(gens(OO(U21))[1])*L(U21)[1]
  W = PrincipalOpenSubset(U[2], [complement_equation(U21), gens(OO(U[2]))[2]])
  rr = L(U[1], W)
  rrr = L(U21, W)
  @test rr == compose(rho, rrr)
  WW = simplify(W)
  @test WW isa Oscar.SimplifiedSpec
  rrWW = L(U[1], WW)
  rrrWW = L(U21, WW)
  @test rrWW == compose(rho, rrrWW)

  M1 = oscar.cotangent_sheaf(X)
  rho = M1(U[1], U21)
  @test rho(M1(U[1])[1]) in M1(U21)
  T = oscar.tangent_sheaf(X)
  rho = T(U[1], U21)
  g = rho(T(U[1])[1]) 
  @test g in T(U21)
  @test element_to_homomorphism(g)(domain(T)(U21)[1]) in codomain(T)(U21)

  simplify!(X)
  CC = coverings(X)[2]
  for U in patches(CC)
    @test !iszero(T(U))
  end
  W = PrincipalOpenSubset(U[1], one(OO(U[1])))
  @test !iszero(T(W))

  HomM1M1 = oscar.HomSheaf(M1, M1)
  rho = HomM1M1(U[1], U21)
  g = HomM1M1(U[1])[1]
  @test rho(g) in HomM1M1(U21)
  @test element_to_homomorphism(rho(g))(domain(HomM1M1)(U21)[1]) in codomain(HomM1M1)(U21)
end

@testset "Pushforward of modules" begin
  IP = projective_space(QQ, ["x", "y", "z"])
  S = ambient_coordinate_ring(IP)
  (x,y,z) = gens(S)
  f = z^2 - x*y
  I = ideal(S, f)
  II = IdealSheaf(IP, I)
  X = covered_scheme(IP)
  inc = Oscar.CoveredClosedEmbedding(X, II)
  C = domain(inc)
  TC = tangent_sheaf(C)
  U = affine_charts(C)
  TC1 = TC(U[1])
  U21 = PrincipalOpenSubset(U[2], gens(OO(C)(U[2]))[1])
  U321 = PrincipalOpenSubset(U[3], gens(OO(C)(U[3]))[1]*gens(OO(U[3]))[2])
  @test !iszero(TC(U321))
  @test !iszero(TC(U21))
  @test TC(U[1], U321)(TC1[1]) == TC(U21, U321)(TC(U[1], U21)(TC1[1]))
  incTC = Oscar.PushforwardSheaf(inc, TC)
  U = affine_charts(X)
  TC1 = incTC(U[1])
  @test TC1 === incTC(U[1])
  @test incTC(U[1]) isa SubQuo
  U21 = PrincipalOpenSubset(U[2], dehomogenize(IP, 1)(x))
  U321 = PrincipalOpenSubset(U[3], dehomogenize(IP, 2)(x*y))
  @test incTC(U[1], U321)(TC1[1]) == incTC(U21, U321)(incTC(U[1], U21)(TC1[1]))
end
