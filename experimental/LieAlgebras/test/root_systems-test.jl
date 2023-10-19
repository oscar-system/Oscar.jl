@testset "LieAlgebras.RootSystem" begin
  R = root_system(:F, 4)
  @test cartan_matrix(R) ==
    matrix(QQ, 4, 4, [2, -1, 0, 0, -1, 2, -2, 0, 0, -1, 2, -1, 0, 0, -1, 2])
  @test number_of_roots(R) == 48
  @test number_of_roots(:F, 4) == number_of_roots(R)
end
