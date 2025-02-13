@testset "Schur polynomials" begin
  S,x = polynomial_ring(ZZ, :x => 1:3)

  # schur_polynomial(R, λ)
  @test schur_polynomial(S, partition([1])) == x[1]
  @test schur_polynomial(S, partition([1,1])) == x[1]*x[2]
  @test schur_polynomial(S, partition([2])) == x[1]^2
  @test schur_polynomial(S, partition([1,1,1])) == x[1]*x[2]*x[3]
  @test schur_polynomial(S, partition([2,1])) == x[1]^2*x[2] + x[1]*x[2]^2
  @test schur_polynomial(S, partition([3])) == x[1]^3
  @test schur_polynomial(S, partition([])) == 1

  # schur_polynomial(R, λ, n)
  @test schur_polynomial(S, partition([1]), 3) == x[1] + x[2] + x[3]
  @test schur_polynomial(S, partition([1,1]), 3) == x[1]*x[2] + x[1]*x[3] + x[2]*x[3]
  @test schur_polynomial(S, partition([2]), 2) == x[1]^2 + x[2]^2 + x[1]*x[2]
  @test schur_polynomial(S, partition([2]), 3) == x[1]^2 + x[2]^2 + x[1]*x[2] + x[1]*x[3] + x[2]*x[3] + x[3]^2
  @test schur_polynomial(S, partition([1,1,1]), 3) == x[1]*x[2]*x[3]
  @test schur_polynomial(S, partition([2,1]), 1) == 0
  @test schur_polynomial(S, partition([2,1]), 3) == 2*x[1]*x[2]*x[3] + x[1]^2*x[2] + x[1]*x[2]^2 + x[1]^2*x[3] + x[1]*x[3]^2 + x[2]^2*x[3] + x[2]*x[3]^2
  @test schur_polynomial(S, partition([3,1]), 2) == x[1]^3*x[2] + x[1]^2*x[2]^2 + x[1]*x[2]^3
  @test schur_polynomial(S, partition([3]), 2) == x[1]^2*x[2] + x[1]*x[2]^2 + x[1]^3 + x[2]^3
  @test schur_polynomial(S, partition([3]), 3) == x[1]*x[2]*x[3] + x[1]^2*x[2] + x[1]*x[2]^2 + x[1]^2*x[3] + x[1]*x[3]^2 + x[2]^2*x[3] + x[2]*x[3]^2 + x[1]^3 + x[2]^3 + x[3]^3

  @test schur_polynomial(S, partition([141]), 1) == x[1]^141 #this calls Cauchy's algorithm
  @test schur_polynomial(S, partition([200]), 1) == x[1]^200 #this calls Cauchy's algorithm

  @test schur_polynomial(S, partition([]), 1) == 1
  @test schur_polynomial(S, partition([]), 0) == 1
  @test schur_polynomial(S, partition([3,2,1]), 0) == 0
  @test_throws ArgumentError schur_polynomial(partition([3,2,1]), -1)

  #schur_polynomial(λ)
  @test schur_polynomial(partition([])) == 1
  @test schur_polynomial(partition([1]))(1) == 1

  #schur_polynomial(λ, n)
  @test schur_polynomial(partition([]), 1) == 1
  @test schur_polynomial(partition([]), 0) == 1
  @test schur_polynomial(partition([3,2,1]), 0) == 0
  @test schur_polynomial(partition([1]), 10)(1,1,1,1,1,1,1,1,1,1) == 10
  @test_throws ArgumentError schur_polynomial(partition([3,2,1]), -1)

  #Two examples from Wikipedia
  schur_polynomial(S, partition([2,1,1]), 3) == x[1]*x[2]*x[3]*(x[1]+x[2]+x[3])
  schur_polynomial(S, partition([2,2]), 3) == x[1]^2*x[2]^2 + x[1]^2*x[3]^2 + x[2]^2*x[3]^2 + x[1]^2*x[2]*x[3] + x[1]*x[2]^2*x[3] + x[1]*x[2]*x[3]^2

  # From issue #3850
  R, _ = polynomial_ring(ZZ, [:x1, :x2, :x3])
  f = Oscar.schur_polynomial_cbf(R, partition([49, 2]))
  g = Oscar.schur_polynomial_combinat(R, partition([49, 2]))
  @test f == g
end
