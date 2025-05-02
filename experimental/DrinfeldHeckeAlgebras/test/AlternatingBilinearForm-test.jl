@testset "DrinfeldHeckeAlgebras.AlternatingBilinearForm" begin
  MS = matrix_space(QQ,2,2)
  m1 = MS([1 2;3 4])
  m2 = MS([0 1;-1 0])
  
  b1 = BilinearForm(m1)
  b2 = alternating_bilinear_form(m2)
  
  @test matrix(b1) == m1
  @test matrix(b2) == m2
  @test_throws ArgumentError alternating_bilinear_form(m1)
  
  v = [QQ(1), QQ(2)]
  w = [QQ(-1), QQ(3)]
  
  @test b1(v,w) == 23
  @test b2(v,w) == 5
end
