@testset "LieAlgebras.Combinatorics" begin
  @testset "permutations" begin
    permutations = Oscar.LieAlgebras.permutations

    for n in 1:6
      @test collect(permutations(n)) == collect(permutations(1:n))
    end

    @test Set(permutations(["1", "2", "3"])) == Set([
      ["1", "2", "3"],
      ["1", "3", "2"],
      ["2", "1", "3"],
      ["2", "3", "1"],
      ["3", "1", "2"],
      ["3", "2", "1"],
    ])
  end

  @testset "permutations_with_sign" begin
    permutations_with_sign = Oscar.LieAlgebras.permutations_with_sign

    for n in 1:6
      @test collect(permutations_with_sign(n)) == collect(permutations_with_sign(1:n))
    end

    @test Set(permutations_with_sign(["1", "2", "3"])) == Set([
      (["1", "2", "3"], 1),
      (["1", "3", "2"], -1),
      (["2", "1", "3"], -1),
      (["2", "3", "1"], 1),
      (["3", "1", "2"], 1),
      (["3", "2", "1"], -1),
    ])
  end
end
