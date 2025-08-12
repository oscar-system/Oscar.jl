@testset "Linear quotients" begin
  G = matrix_group(matrix(QQ, 1, 1, [2]))
  @test_throws ArgumentError linear_quotient(G)

  K, a = cyclotomic_field(4, "a")
  r = matrix(K, 2, 2, [-1, 0, 0, 1])
  s = matrix(K, 2, 2, [a, 0, 0, inv(a)])
  G = matrix_group(r, s)
  L = linear_quotient(G)

  @test group(L) === G
  @test base_ring(L) === K

  A, GtoA = class_group(L)
  @test is_snf(A)
  @test A.snf == ZZRingElem[2]
  @test is_zero(GtoA(G(r)))
  @test GtoA(G(s)) == A([1])

  K, a = cyclotomic_field(6, "a")
  g = diagonal_matrix([a^2, a^2, a^2])
  G = matrix_group(g)
  @test Oscar.age(G(g), (a, 6)) == 1
  @test Oscar.age(G(g), (a^2, 3)) == 1
  @test Oscar.age(G(g), (a^5, 6)) == 2
  @test Oscar.representatives_of_junior_elements(G, (a^2, 3)) == elem_type(G)[G(g)]
  L = linear_quotient(G)
  @test has_canonical_singularities(L)
  @test !has_terminal_singularities(L)

  g = diagonal_matrix([QQ(-1), QQ(-1), QQ(-1), QQ(-1)])
  G = matrix_group(g)
  @test Oscar.age(G(g), (QQ(-1), 2)) == 2
  @test isempty(Oscar.representatives_of_junior_elements(G, (QQ(-1), 2)))
  L = linear_quotient(G)
  @test has_canonical_singularities(L)
  @test has_terminal_singularities(L)

  g = diagonal_matrix([QQ(-1)])
  G = matrix_group(g)
  @test Oscar.age(G(g), (QQ(-1), 2)) == QQFieldElem(1, 2)
  L = linear_quotient(G)
  @test !has_canonical_singularities(L)
  @test !has_terminal_singularities(L)

  K, a = cyclotomic_field(12, "a")
  g = matrix(K, 2, 2, [0, a^4, a^4, 0])
  G = matrix_group(g)
  R, x = polynomial_ring(K, 2)
  val = Oscar.monomial_valuation(R, G(g), (a, 12))
  @test val(x[1]) == 2
  @test val(x[2]) == 2
  @test val(x[1] + x[2]) == 2
  @test val(x[1] - x[2]) == 5
  val = Oscar.monomial_valuation(R, G(g), (a^2, 6))
  @test val(x[1]) == 2
  @test val(x[2]) == 2
  @test val(x[1] + x[2]) == 2
  @test val(x[1] - x[2]) == 5
  val = Oscar.monomial_valuation(R, G(g), (a^-1, 12))
  @test val(x[1]) == 1
  @test val(x[2]) == 1
  @test val(x[1] + x[2]) == 4
  @test val(x[1] - x[2]) == 1

  # An example containing pseudo-reflections
  K, a = cyclotomic_field(6, "a")
  g1 = matrix(K, 3, 3, [-1 0 0; 0 1 0; 0 0 1])
  g2 = matrix(K, 3, 3, [0 0 1; 1 0 0; 0 1 0])
  G = matrix_group(g1, g2)
  L = linear_quotient(G)
  A, GtoA = class_group(L)
  @test is_snf(A)
  @test A.snf == ZZRingElem[3]
  @test is_zero(GtoA(G(g1)))
  @test GtoA(G(g2)) == A([1])
end
