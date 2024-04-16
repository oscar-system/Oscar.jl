@testset "Cox rings of linear quotients" begin
  C2 = matrix_group(matrix(QQ, 2, 2, [-1, 0, 0, -1]))
  L = linear_quotient(C2)
  R, RtoS = cox_ring(L)
  @test is_zero(modulus(R))
  @test ngens(R) == 2
  A = grading_group(R)
  @test A.snf == ZZRingElem[2]
  @test degree(gen(R, 1)) == A([1])
  @test degree(gen(R, 2)) == A([1])

  K, a = cyclotomic_field(4, "a")
  D4 = matrix_group(matrix(K, 2, 2, [a, 0, 0, inv(a)]), matrix(K, 2, 2, [0, -1, 1, 0]))
  L = linear_quotient(D4)
  R, RtoS = cox_ring(L)
  A, GtoA = class_group(L)
  @test A === grading_group(R)
  for i in 1:ngens(R)
    @test Oscar.ab_g_degree(GtoA, RtoS(gen(R, i)), Oscar.fixed_root_of_unity(L)) ==
      degree(gen(R, i))
  end

  K, a = cyclotomic_field(60, "a")
  b = a^15
  M1 = matrix(K, 2, 2, [b, 0, 0, inv(b)])
  c = a^6
  r = c + inv(c)
  s = b * (c^3 + inv(c)^3)
  M2 = K(1//2) * matrix(K, 2, 2, [r + s, K(1), K(-1), r - s])
  I = matrix_group([M1, M2])
  L = linear_quotient(I)
  R, RtoS = cox_ring(L; algo_rels=:linear_algebra)
  @test is_trivial(grading_group(R))
  A, GtoA = class_group(L)
  @test A === grading_group(R)
  for i in 1:ngens(R)
    @test Oscar.ab_g_degree(GtoA, RtoS(gen(R, i)), Oscar.fixed_root_of_unity(L)) ==
      degree(gen(R, i))
  end

  K, a = cyclotomic_field(12, "a")
  g1 = matrix(K, 4, 4, [1 0 0 0; 0 a^4 0 0; 0 0 1 0; 0 0 0 a^-4])
  g2 =
    1//3 * matrix(
      K,
      4,
      4,
      [
        2 * a^4+1 a^4-1 0 0
        2 * a^4-2 a^4+2 0 0
        0 0 -2 * a^4-1 -2 * a^4-4
        0 0 -a^4-2 -a^4+1
      ],
    )
  G = matrix_group(g1, g2) # symplectic reflection representation of reflection group G_4
  L = linear_quotient(G)
  R, RtoS = cox_ring(L; algo_rels=:linear_algebra)
  @test ngens(R) == 18
  @test grading_group(R).snf == ZZRingElem[3]
  A, GtoA = class_group(L)
  @test A === grading_group(R)
  for i in 1:ngens(R)
    @test Oscar.ab_g_degree(GtoA, RtoS(gen(R, i)), Oscar.fixed_root_of_unity(L)) ==
      degree(gen(R, i))
  end

  # An example containing reflections
  K, a = cyclotomic_field(6, "a")
  g1 = matrix(K, 3, 3, [-1 0 0; 0 1 0; 0 0 1])
  g2 = matrix(K, 3, 3, [0 0 1; 1 0 0; 0 1 0])
  G = matrix_group(g1, g2)
  L = linear_quotient(G)
  R, RtoS = cox_ring(L)
  @test ngens(R) == 3
  @test grading_group(R).snf == ZZRingElem[3]
  A, GtoA = class_group(L)
  @test A === grading_group(R)
  for i in 1:ngens(R)
    @test Oscar.ab_g_degree(GtoA, RtoS(gen(R, i)), Oscar.fixed_root_of_unity(L)) ==
      degree(gen(R, i))
  end
  @test is_zero(modulus(R))
end

@testset "Cox rings of QQ-factorial terminalizations" begin
  K, a = cyclotomic_field(4)
  G = matrix_group(
    matrix(K, 4, 4, [0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0]),
    matrix(K, 4, 4, [a, 0, 0, 0, 0, inv(a), 0, 0, 0, 0, inv(a), 0, 0, 0, 0, a]),
  ) # symplectic reflection representation of dihedral group of order 8
  L = linear_quotient(G)
  C, CtoRt = Oscar.cox_ring_of_qq_factorial_terminalization(L)
  # Can't really do tests of correctness, so this is more a test of "does the
  # function still exists"
  @test is_isomorphic(grading_group(C), abelian_group([0, 0]))
  for i in 1:ngens(C)
    f = CtoRt(gen(C, i))
    exps = collect(AbstractAlgebra.exponent_vectors(f))
    @test length(exps) == 1
    @test degree(gen(C, i)) == grading_group(C)(exps[1])
  end
end
