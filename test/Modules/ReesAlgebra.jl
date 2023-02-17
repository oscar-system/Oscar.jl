@testset "Rees algebras" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  F = FreeMod(R, 3)
  M, inc = sub(F, [a*g for a in gens(R) for g in gens(F)])
  Q, p = quo(F, [a*F[1] for a in gens(R)])

  S1 = oscar.rees_algebra(F)
  @test !iszero(one(S1))
  S2 = oscar.rees_algebra(M)
  @test !iszero(one(S2))
  S3 = oscar.rees_algebra(Q)
  @test !iszero(one(S3))

  A = R[x y; x-1 z]
  f = det(A)
  I = ideal(R, f)
  Q, _ = quo(R, I)
  AQ = change_base_ring(Q, A)
  FQ = FreeMod(Q, 2)
  M, _ = quo(FQ, change_base_ring(Q, A))
  M_double_dual, psi = oscar.double_dual(M)

  @test is_isomorphism(psi)
  f_A = hom(FQ, FQ, AQ)
  dual_f_A = dual(f_A)
  @test codomain(dual_f_A) isa FreeMod
  g = oscar._versal_morphism_to_free_module(M)

  RM = oscar.rees_algebra(g)
  @test !iszero(one(RM))
end
