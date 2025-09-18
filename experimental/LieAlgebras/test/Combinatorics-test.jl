@testset "LieAlgebras.Combinatorics" begin
  @testset "multicombinations" begin
    multicombinations = Oscar.LieAlgebras.multicombinations

    @test collect(multicombinations(0, 0)) == [Int[]]
    @test collect(multicombinations(0, 1)) == []

    for n in 1:8, k in 1:n
      @test collect(multicombinations(n, k)) == collect(multicombinations(1:n, k))
    end

    @test collect(multicombinations(["1", "2", "3", "4"], 0)) == [String[]]
    @test collect(multicombinations(["1", "2", "3", "4"], 1)) ==
      [["1"], ["2"], ["3"], ["4"]]
    @test collect(multicombinations(["1", "2", "3", "4"], 2)) == [
      ["1", "1"],
      ["1", "2"],
      ["1", "3"],
      ["1", "4"],
      ["2", "2"],
      ["2", "3"],
      ["2", "4"],
      ["3", "3"],
      ["3", "4"],
      ["4", "4"],
    ]
    @test collect(multicombinations(["1", "2", "3", "4"], 3)) == [
      ["1", "1", "1"],
      ["1", "1", "2"],
      ["1", "1", "3"],
      ["1", "1", "4"],
      ["1", "2", "2"],
      ["1", "2", "3"],
      ["1", "2", "4"],
      ["1", "3", "3"],
      ["1", "3", "4"],
      ["1", "4", "4"],
      ["2", "2", "2"],
      ["2", "2", "3"],
      ["2", "2", "4"],
      ["2", "3", "3"],
      ["2", "3", "4"],
      ["2", "4", "4"],
      ["3", "3", "3"],
      ["3", "3", "4"],
      ["3", "4", "4"],
      ["4", "4", "4"],
    ]
  end

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
