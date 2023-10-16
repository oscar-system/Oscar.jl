@testset "affine algebraic sets" begin
  P, (x,y) = polynomial_ring(QQ, [:x, :y])
  p = x^2 - y^2
  X =  algebraic_set(p)
  C = irreducible_components(X)
  @test length(C) == 2
  Y = algebraic_set(p^2)
  @test X == Y

  X1 =  algebraic_set(x^2 + y^2) # irreducible but not geometrically irreducible.
  @test length(irreducible_components(X1)) == 1
  C1 = geometric_irreducible_components(X1)
  @test length(C1) == 1
  @test C1[1][3]==2
end
