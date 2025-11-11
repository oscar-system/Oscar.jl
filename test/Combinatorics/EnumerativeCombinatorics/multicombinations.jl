@testset "multicombinations" begin
  let
    C = multicombinations(3, 2)
    @test length(C) == 6
    @test collect(C) == [[1, 1],
                         [1, 2],
                         [1, 3],
                         [2, 2],
                         [2, 3],
                         [3, 3]]
    @test eltype(C) === Combination{Int}
    C = multicombinations(100, 3)
    @test length(C) == binomial(102, 3)

    C = multicombinations(Int16(2), Int16(2))
    @test length(C) == 3
    @test collect(C) == [[1, 1],
                         [1, 2],
                         [2, 2]]
    @test eltype(C) === Combination{Int16}

    C = multicombinations('a':'c', 3)
    @test length(C) == 10
    @test collect(C) == [['a', 'a', 'a'],
                         ['a', 'a', 'b'],
                         ['a', 'a', 'c'],
                         ['a', 'b', 'b'],
                         ['a', 'b', 'c'],
                         ['a', 'c', 'c'],
                         ['b', 'b', 'b'],
                         ['b', 'b', 'c'],
                         ['b', 'c', 'c'],
                         ['c', 'c', 'c']]
    @test eltype(C) === Combination{Char}
  end

  for n in 1:5, k in 1:n
    @test collect(multicombinations(n, k)) == collect(multicombinations(1:n, k))
  end

  @test collect(multicombinations(["1", "2", "3", "4"], 0)) == [String[]]
  @test collect(multicombinations(["1", "2", "3", "4"], 1)) == [["1"], ["2"], ["3"], ["4"]]
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

  v = collect(multicombinations(5, 3))
  @test issorted(v)
  @test allunique(v)

  @test collect(multicombinations(0, 0)) == [Int[]]
  @test collect(multicombinations(0, 1)) == []

  @testset "inplace iteration" begin
    Ci = multicombinations(5,3, inplace=true)
    Cf = multicombinations(5,3)
    @test all(splat(==), zip(Ci, Cf))
  end

end
