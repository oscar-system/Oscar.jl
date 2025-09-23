@testset "Finite fields, discrete logarithm" begin
  F, z = finite_field(2, 4)
  @test is_primitive(z^2)
  @test !is_primitive(z^3)
  @test is_primitive(Oscar.DiscLog.generator(F))
  @test [disc_log(z, z^i) for i in 1:14] == 1:14
  @test disc_log(z, z^15) == 0
  @test disc_log(z^2, z) == 8
  @test_throws "disc_log failed" disc_log(z, zero(z))
end
