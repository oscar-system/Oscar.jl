@testset "HasseSchmidt Derivatives MPolyRingElem" begin
  R, (x, y) = polynomial_ring(ZZ, ["x", "y"]);
  
  @test [x^3, 3x^2, 3x, R(1)] == hasse_derivatives(R(x^3))
  @test [5x^2 + 3y^5, 15y^4, 30y^3, 30y^2, 15y, R(3), 10x, R(5)] == hasse_derivatives(R(5*x^2 + 3*y^5))
  @test [[3*x^2, 6*x, R(3)], [6*y^3, 18*y^2, 18*y, R(6)]] == hasse_derivatives([R(3*x^2), R(6*y^3)])
end

@testset "HasseSchmidt Derivatives MPolyQuoRingElem" begin
  


end

@testset "HasseSchmidt Derivatives Oscar.MPolyLocRingElem" begin
  

  
end

@testset "HasseSchmidt Derivatives Oscar.MPolyQuoLocRingElem" begin
  

  
end