@testset "Experimental.DiscreteLog" begin

  F = GF(3,4);
  a = gen(F)^21;
  @test Oscar.DiscreteLog.disc_log(gen(F), a) == 21
  @test_throws "disc_log failed" Oscar.DiscreteLog.disc_log(one(F), a)

  F2 = GF(3,5);
  @test_throws AssertionError Oscar.DiscreteLog.disc_log(gen(F), gen(F2))
end
