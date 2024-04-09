@testset "weak compositions" begin
  @test length(weak_compositions(3, 2)) == 4
  @test @inferred collect(weak_compositions(3, 2)) == map(weak_composition, [[3, 0], [2, 1], [1, 2], [0, 3]])
  @test @inferred collect(weak_compositions(Int8(3), 2)) == map(weak_composition, Vector{Int8}[[3, 0], [2, 1], [1, 2], [0, 3]])
  @test isempty(weak_compositions(1, 0))
  @test @inferred collect(weak_compositions(1, 0)) == Oscar.WeakComposition{Int}[]
  @test length(weak_compositions(0, 0)) == 1
  @test @inferred collect(weak_compositions(0, 0)) == [weak_composition(Int[])]
  @test length(weak_compositions(0, 3)) == 1
  @test @inferred collect(weak_compositions(0, 3)) == [weak_composition(Int[0, 0, 0])]
end
