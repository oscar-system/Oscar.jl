@testset "koszul complexes" begin
  S, x = graded_polynomial_ring(QQ, [:x, :y, :z])

  r = length(x)
  d = 5
  F = (is_graded(S) ? graded_free_module(S, degree.(x)) : FreeMod(S, r))
  v = sum(u^d*e for (u, e) in zip(x, gens(F)); init = zero(F))

  K = koszul_complex(Oscar.KoszulComplex, v)
  KK = koszul_complex(v)

  for i in 0:4
    @test matrix(map(K, i)) == matrix(map(KK, i))
  end
  
  # the ungraded case
  S, x = polynomial_ring(QQ, [:x, :y, :z])

  r = length(x)
  d = 5
  F = (is_graded(S) ? graded_free_module(S, degree.(x)) : FreeMod(S, r))
  v = sum(u^d*e for (u, e) in zip(x, gens(F)); init = zero(F))

  K = koszul_complex(Oscar.KoszulComplex, [u^d for u in gens(S)]...)
  KK = koszul_complex(v)

  for i in 0:4
    @test matrix(map(K, i)) == matrix(map(KK, i))
  end
end

@testset "induced morphisms on Koszul complexes" begin
  R, (x, y, z) = QQ[:x, :y, :z]
  F = FreeMod(R, 3)
  v = x^3*F[1] + y^2*F[2] + z^4*F[3]

  G  = FreeMod(R, 2)
  w = (x^3 + 3*y^2)*G[1] + z^4*G[2]

  phi = hom(F, G, [G[1], 3*G[1], G[2]])

  KF = koszul_complex(Oscar.KoszulComplex, v)
  KG = koszul_complex(Oscar.KoszulComplex, w)

  f = koszul_complex(phi, KF, KG)
  @test f[1] isa FreeModuleHom
  @test f[2] isa FreeModuleHom
  @test f[3] isa FreeModuleHom
  @test compose(f[3], map(KG, 2)) == compose(map(KF, 3), f[2])
  @test compose(f[2], map(KG, 1)) == compose(map(KF, 2), f[1])
end

@testset "homogeneous Koszul complexes and their induced maps" begin
  S, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])

  a = [x, y, z]
  b = [x^2, y^3, y^5+z^5, x+y]

  KA = Oscar.HomogKoszulComplex(S, a)
  KB = Oscar.HomogKoszulComplex(S, b)

  phi = Oscar.InducedKoszulMorphism(KB, KA)

  G = grading_group(S)
  S1 = graded_free_module(S, [zero(G)])
  KAD = hom(KA, S1)
  KBD = hom(KB, S1)

  phi_star = hom(phi, Oscar.ZeroDimensionalComplex(S1))
  phi_star[0]
  @test_throws ErrorException phi_star[1]
  phi_star[-1]
  phi_star[-2]
  phi_star[-3]

  @test compose(phi_star[-1], map(codomain(phi_star), -1)) == compose(map(domain(phi_star), -1), phi_star[-2])
end

