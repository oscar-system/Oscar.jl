@testset "compositions" begin
  compositions = Oscar.compositions
  @test isempty(compositions(1, 0))
  @test length(compositions(0, 0)) == 1
  @test isempty(compositions(2, 3))
  @test [1, 2] in compositions(3, 2)
  @test [2, 1] in compositions(3, 2)
  @test [1, 3] in compositions(4, 2)
  @test [2, 2] in compositions(4, 2)
  @test [3, 1] in compositions(4, 2)
end
