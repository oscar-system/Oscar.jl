@testset "constant sheaf of integers on AffineScheme" begin
  R, x = QQ[:x, :y, :z]
  X = spec(R)
  Xcov = covered_scheme(X)

  mutable struct ZZConstantSheaf{SpaceType, OpenType, OutputType,
                                 RestrictionType
                                } <: AbsPreSheaf{
                                                 SpaceType, OpenType,
                                                 OutputType, RestrictionType
                                                }
    und::PreSheafOnScheme
    function ZZConstantSheaf(und::PreSheafOnScheme{A, B, C, D}) where {A, B, C, D}
      return new{A, B, C, D}(und)
    end
  end

  Oscar.underlying_presheaf(F::ZZConstantSheaf) = F.und

  function Oscar.produce_object(F::ZZConstantSheaf, U::AbsAffineScheme)
    return ZZ
  end

  function Oscar.produce_restriction_map(F::ZZConstantSheaf, V::AbsAffineScheme, U::AbsAffineScheme)
    return identity_map(ZZ)
  end

  und = PreSheafOnScheme(X, 
                         OpenType=AbsAffineScheme, OutputType=typeof(ZZ),
                         RestrictionType=typeof(identity_map(ZZ)),
                         is_open_func=is_open_embedding
                        )
  const_sheaf_ZZ = ZZConstantSheaf(und)

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
  P = projective_space(QQ, 2)
  PC = covered_scheme(P)
  F = StructureSheafOfRings(PC)
  U = patches(default_covering(PC))
  R0 = F(U[1])
  yx, zx = gens(R0)
  W = PrincipalOpenSubset(U[1], yx)
  rho = F(U[2], W)
  @test rho === F(U[2], W)
  R1 = F(U[2])
  xy, zy = gens(R1)
  @test rho(xy) == inv(OO(W)(yx))
  @test rho(zy) == OO(W)(zx)*inv(OO(W)(yx))
  rho2 = F(U[1], W)
  @test rho2.(gens(OO(U[1]))) == OO(W).(gens(OO(U[1])))

  B = PrincipalOpenSubset(U[2], gens(OO(U[2]))[1])
  eta = F(U[1], B)
  @test eta === F(U[1], B)
  tmp1 = F(W, B)
  tmp2 = F(B, W)
  @test tmp1(gens(OO(W))[1])== inv(gens(OO(B))[1])
  @test tmp2(gens(OO(B))[1])== inv(gens(OO(W))[1])
  h = compose(tmp1, tmp2)
  @test h(gens(OO(W))[1]) == gens(OO(W))[1]

  O = AffineSchemeOpenSubscheme(U[1], lifted_numerator.([yx*(yx-1), yx*(zx-1)]))
  @test F(U[1], O)(OO(U[1])[1]) == OO(O)(OO(U[1])[1])
  @test F(U[2], O)(OO(U[2])[1]) == inv(OO(O)(OO(U[1])[1]))

  # To test the other gluings, we repeat the above with a new gluing.
  P = projective_space(QQ, 2)
  PC = covered_scheme(P)
  cov = default_covering(PC)
  U = patches(cov)
  for i in 1:length(U)-1
    for j in i+1:length(U)
      G0 = cov[U[i], U[j]]
      G1 = Gluing(G0)
      add_gluing!(cov, G1)
    end
  end

  F = StructureSheafOfRings(PC)
  R0 = F(U[1])
  yx, zx = gens(R0)
  W = PrincipalOpenSubset(U[1], yx)
  @test F(U[1], W)(gen(R0, 1)) == OO(W)(gen(R0, 1))
  rho = F(U[2], W)
  @test rho === F(U[2], W)
  R1 = F(U[2])
  xy, zy = gens(R1)
  @test rho(xy) == inv(OO(W)(yx))
  @test rho(zy) == OO(W)(zx)*inv(OO(W)(yx))
  rho2 = F(U[1], W)
  @test rho2.(gens(OO(U[1]))) == OO(W).(gens(OO(U[1])))

  B = PrincipalOpenSubset(U[2], gens(OO(U[2]))[1])
  eta = F(U[1], B)
  @test eta === F(U[1], B)
  tmp1 = F(W, B)
  tmp2 = F(B, W)
  @test tmp1(gens(OO(W))[1])== inv(gens(OO(B))[1])
  @test tmp2(gens(OO(B))[1])== inv(gens(OO(W))[1])
  h = compose(tmp1, tmp2)
  @test h(gens(OO(W))[1]) == gens(OO(W))[1]

  O = AffineSchemeOpenSubscheme(U[1], lifted_numerator.([yx*(yx-1), yx*(zx-1)]))
  @test F(U[1], O)(OO(U[1])[1]) == OO(O)(OO(U[1])[1])
  @test F(U[2], O)(OO(U[2])[1]) == inv(OO(O)(OO(U[1])[1]))

  G = cov[U[1], U[2]]
  f, g = gluing_morphisms(G)
  O2 = preimage(g, O)
  @test F(O2, O)(gens(OO(O2))[1]) == inv(OO(O)(OO(U[1])[1]))
  @test F(O, O2)(inv(OO(O)(OO(U[1])[1]))) == gens(OO(O2))[1]
  O3 = AffineSchemeOpenSubscheme(U[1], lifted_numerator.([yx*(yx-1)^3, yx*(zx-1)]))
  @test F(O, O3)(one(OO(O))) == one(OO(O3))
  @test F(O2, O3)(one(OO(O2))) == one(OO(O3))
end
