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

@testset "weak compositions" begin
  weak_compositions = Oscar.weak_compositions
  @test isempty(weak_compositions(1, 0))
  @test length(weak_compositions(0, 0)) == 1
  @test [0, 3] in weak_compositions(3, 2)
  @test [3, 0] in weak_compositions(3, 2) 
  @test length(weak_compositions(1, 3)) == 3  
end
