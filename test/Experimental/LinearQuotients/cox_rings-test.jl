@testset "Cox rings of linear quotients" begin
  C2 = matrix_group(matrix(QQ, 2, 2, [ -1, 0, 0, -1 ]))
  L = linear_quotient(C2)
  R, RtoS = cox_ring(L)
  @test is_zero(modulus(R))
  @test ngens(R) == 2
  A = grading_group(R)
  @test A.snf == ZZRingElem[ 2 ]
  @test degree(gen(R, 1)) == A([ 1 ])
  @test degree(gen(R, 2)) == A([ 1 ])

  K, a = cyclotomic_field(4, "a")
  D4 = matrix_group(matrix(K, 2, 2, [ a, 0, 0, inv(a) ]), matrix(K, 2, 2, [ 0, -1, 1, 0 ]))
  L = linear_quotient(D4)
  R, RtoS = cox_ring(L)
  A, GtoA = class_group(L)
  @test A === grading_group(R)
  for i = 1:ngens(R)
    @test Oscar.ab_g_degree(GtoA, RtoS(gen(R, i)), Oscar.fixed_root_of_unity(L)) == degree(gen(R, i))
  end

  K, a = CyclotomicField(60, "a")
  b = a^15
  M1 = matrix(K, 2, 2, [ b, 0, 0, inv(b) ])
  c = a^6
  r = c + inv(c)
  s = b*(c^3 + inv(c)^3)
  M2 = K(1//2)*matrix(K, 2, 2, [ r + s, K(1), K(-1), r - s ])
  I = matrix_group([ M1, M2 ])
  L = linear_quotient(I)
  R, RtoS = cox_ring(L, algo_rels = :linear_algebra)
  @test is_trivial(grading_group(R))
  A, GtoA = class_group(L)
  @test A === grading_group(R)
  for i = 1:ngens(R)
    @test Oscar.ab_g_degree(GtoA, RtoS(gen(R, i)), Oscar.fixed_root_of_unity(L)) == degree(gen(R, i))
  end
end
