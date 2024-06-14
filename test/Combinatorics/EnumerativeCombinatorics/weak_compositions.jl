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

  l = @inferred collect(weak_compositions(5, 3))
  @test length(l) == 21
  @test all(x -> sum(x) == 5, l)
  b = unique!(l)
  @test length(b) == length(l)

  l = @inferred collect(weak_compositions(5, 0))
  @test l == Vector{Oscar.WeakComposition{Int}}()
  l = @inferred collect(weak_compositions(0, 7))
  @test l == [weak_composition([0 for i in 1:7])]
  l = @inferred collect(weak_compositions(1, 7))
  @test length(l) == length(unique!(l)) == 7

  l = @inferred collect(weak_compositions(7, 1))
  @test l == [weak_composition([7])]

  @test number_of_weak_compositions(1, -2) == ZZ(0)
  @test number_of_weak_compositions(-1, 0) == ZZ(0)
end
