@testset "toric schemes" begin
  S = hirzebruch_surface(3)
  X = ToricCoveredScheme(S)
  @test issmooth(X)

  IP1 = projective_space(NormalToricVariety, 1)
  IP1xIP1 = IP1*IP1
  Y = ToricCoveredScheme(IP1xIP1)
  @test length(values(glueings(default_covering(Y)))) == 16
end
