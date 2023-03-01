@testset "Linear quotients" begin
  G = matrix_group(matrix(QQ, 1, 1, [ 2 ] ))
  @test_throws AssertionError linear_quotient(G)

  K, a = cyclotomic_field(4, "a")
  r = matrix(K, 2, 2, [ -1, 0, 0, 1 ])
  s = matrix(K, 2, 2, [ a, 0, 0, inv(a) ])
  G = matrix_group(r, s)
  L = linear_quotient(G)

  @test group(L) === G
  @test base_ring(L) === K

  A, GtoA = class_group(L)
  @test is_snf(A)
  @test A.snf == ZZRingElem[ 2 ]
  @test is_zero(GtoA(G(r)))
  @test GtoA(G(s)) == A([ 1 ])
end
