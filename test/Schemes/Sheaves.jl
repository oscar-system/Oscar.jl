@testset "constant sheaf of integers on Spec" begin
  R, x = QQ["x", "y", "z"]
  X = Spec(R)

  is_open_func(U::AbsSpec, V::AbsSpec) = is_open_embedding(U, V)

  function production_func(U::AbsSpec)
    return ZZ
  end

  function restriction_func(V::AbsSpec, U::AbsSpec)
    return identity_map(ZZ)
  end

  const_sheaf_ZZ= SheafOnScheme(X, production_func, restriction_func,
                              OpenType=AbsSpec, OutputType=typeof(ZZ), 
                              RestrictionType=typeof(identity_map(ZZ)),
                              is_open_func=is_open_func
                             )

  U = hypersurface_complement(X, [x[1]])
  G = const_sheaf_ZZ(U)
  @test G === const_sheaf_ZZ(U)
  V = hypersurface_complement(U, [x[2]])
  @test G === const_sheaf_ZZ(V)
  rho = const_sheaf_ZZ(U, V)
  @test rho==identity_map(ZZ)
  @test domain(rho) === const_sheaf_ZZ(U)
  @test codomain(rho) === const_sheaf_ZZ(V)
end

@testset "structure sheaf on covered schemes" begin
  P = projective_space(QQ, 1)
  PC = covered_scheme(P)
  F = RingOfRegularFunctions(PC)
  U = patches(default_covering(PC))
  F(U[1])
  W = PrincipalOpenSubset(U[1], gens(OO(U[1]))[1])
  rho = F(U[2], W)
  @test rho === F(U[2], W)
  B = PrincipalOpenSubset(U[2], gens(OO(U[2]))[1])
  eta = F(U[1], B)
  @test eta === F(U[1], B)
  tmp1 = F(W, B)
  tmp2 = F(B, W)
  @test tmp1(gens(OO(W))[1])== inv(gens(OO(B))[1])
  @test tmp2(gens(OO(B))[1])== inv(gens(OO(W))[1])
  h = compose(tmp1, tmp2)
  @test h(gens(OO(W))[1]) == gens(OO(W))[1]
end
