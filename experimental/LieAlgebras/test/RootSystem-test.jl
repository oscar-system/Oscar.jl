@testset "root_system(cartan_matrix::ZZMatrix)" begin
  R = root_system(:F, 4)
  @test num_positive_roots(R) == 24
  @test num_roots(R) == 48

  R = root_system(:G, 2)
  @test num_positive_roots(R) == 6
  @test num_roots(R) == 12

  @test coefficients(positive_root(R, 1)) == QQ[1 0]
  @test coefficients(positive_root(R, 2)) == QQ[0 1]
end
