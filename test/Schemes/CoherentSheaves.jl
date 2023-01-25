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

@testset "pullbacks of modules" begin
  # Note: The PullbackSheafs are not yet fully functional!!!
  IP = projective_space(QQ, 2)
  S = ambient_coordinate_ring(IP)
  X = covered_scheme(IP)
  (x,y,z) = gens(S)
  f = x^3 + y^3 + z^3
  I = ideal(S, f)
  II = IdealSheaf(IP, I)
  inc = oscar.CoveredClosedEmbedding(X, II)
  C = domain(inc)
  L = twisting_sheaf(IP, 2)
  LC = oscar.PullbackSheaf(inc, L)
  U = affine_charts(C)
  @test LC(U[1]) isa FreeMod
  V = PrincipalOpenSubset(U[1], one(OO(U[1])))
  @test LC(V) isa FreeMod
  rho = LC(U[1], V)
  @test rho isa ModuleFPHom
  VV = simplify(V)
  rho = LC(U[1], VV)
  @test rho isa ModuleFPHom

  # Test pullbacks along blowup maps.
  # This is particularly simple, because the underlying CoveringMorphism 
  # of the projection map is given on the default_covering of its domain.
  J = ideal(S, [x-z, y])
  JJ = IdealSheaf(IP, J)
  JJC = pullback(inc, JJ)
  IP_Bl_C = blow_up(JJC)
  Bl_C = covered_scheme(IP_Bl_C)
  p = covered_projection_to_base(IP_Bl_C)
  p_star = pullback(p)
  p_star_LC = p_star(LC)

  for U in affine_charts(Bl_C)
    @test p_star_LC(U) isa FreeMod
  end

  U = first(affine_charts(Bl_C))
  V = simplify(U)
  @test p_star_LC(V) isa FreeMod
  @test p_star_LC(U, V) isa ModuleFPHom
end

@testset "projectivization of vector bundles" begin
  IP = projective_space(QQ, ["x", "y", "z", "w"])
  S = ambient_coordinate_ring(IP)
  (x,y,z,w) = gens(S)
  f = x^4 + y^4 + z^4 + w^4
  IPC = subscheme(IP, f)
  C = covered_scheme(IPC)
  TC = tangent_sheaf(C)
  IPTC = oscar.projectivization(TC)
  Y = covered_scheme(IPTC)
  @test is_smooth(Y)

  L = twisting_sheaf(IPC, 2)

  PL = oscar.projectivization(L)
  p = covered_projection_to_base(PL)
  @test patches(codomain(covering_morphism(p))) == patches(default_covering(C))
end
