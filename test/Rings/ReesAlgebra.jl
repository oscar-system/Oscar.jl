@testset "Rees algebras" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  F = FreeMod(R, 3)
  M, inc = sub(F, [a*g for a in gens(R) for g in gens(F)])
  Q, p = quo(F, [a*F[1] for a in gens(R)])

  S1 = Oscar.rees_algebra(F)
  @test !iszero(one(S1))
  S2 = Oscar.rees_algebra(M)
  @test !iszero(one(S2))
  S3 = Oscar.rees_algebra(Q)
  @test !iszero(one(S3))

  A = R[x y; x-1 z]
  f = det(A)
  I = ideal(R, f)
  Q, _ = quo(R, I)
  AQ = change_base_ring(Q, A)
  FQ = FreeMod(Q, 2)
  M, _ = quo(FQ, change_base_ring(Q, A))
  RM = Oscar.rees_algebra(M)
  @test is_projective(M)[1]
  _, _, sigma = is_projective(M)
  for f in gens(modulus(RM))
    for v in gens(M)
      c = coordinates(sigma(v))
      # convert SRow to Vector
      w = [c[i] for i in 1:ngens(RM)]
      @test iszero(evaluate(f, w))
    end
  end
  @test RM isa MPolyQuoRing
  M_dual, ev = Oscar.dual(M)
  RM_dual = Oscar.rees_algebra(M_dual)
  @test is_projective(M_dual)[1]
  _, _, sigma = is_projective(M_dual)
  for f in gens(modulus(RM_dual))
    for v in gens(M_dual)
      c = coordinates(sigma(v))
      # convert SRow to Vector
      w = [c[i] for i in 1:ngens(RM_dual)]
      @test iszero(evaluate(f, w))
    end
  end

  @test RM_dual isa MPolyQuoRing
  M_double_dual, psi = Oscar.double_dual(M)

  @test is_isomorphism(psi)
  f_A = hom(FQ, FQ, AQ)
  dual_f_A = dual(f_A)
  @test codomain(dual_f_A) isa FreeMod
  g = Oscar._versal_morphism_to_free_module(M)

  Rg = Oscar.rees_algebra(g)
  @test !iszero(one(Rg))
  # Both algebras should be the same despite being created differently.
  @test all(x->iszero(evaluate(x, gens(RM))), gens(modulus(Rg)))
  @test all(x->iszero(evaluate(x, gens(Rg))), gens(modulus(RM)))
end
