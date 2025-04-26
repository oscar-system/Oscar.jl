@testset "combinations" begin
  let
    C = combinations(4, 3)
    @test length(C) == 4
    @test collect(C) == [[1, 2, 3],
                         [1, 2, 4],
                         [1, 3, 4],
                         [2, 3, 4]]
    @test eltype(C) === Combination{Int}
    C = combinations(100, 3)
    @test length(C) == binomial(100, 3)
    C = combinations(4, 5)
    @test length(C) == 0

    C = combinations('a':'d', 3)
    @test length(C) == 4
    @test collect(C) == [['a', 'b', 'c'],
                         ['a', 'b', 'd'],
                         ['a', 'c', 'd'],
                         ['b', 'c', 'd']]
    @test eltype(C) === Combination{Char}
    C = combinations('a':'d', 10)
    @test length(C) == 0
  end

  for n in 1:8, k in 1:n
    @test collect(combinations(n, k)) == collect(combinations(1:n, k))
  end

  @test collect(combinations(["1", "2", "3", "4"], 0)) == [String[]]
  @test collect(combinations(["1", "2", "3", "4"], 1)) == [["1"], ["2"], ["3"], ["4"]]
  @test collect(combinations(["1", "2", "3", "4"], 2)) ==
  [["1", "2"], ["1", "3"], ["1", "4"], ["2", "3"], ["2", "4"], ["3", "4"]]
  @test collect(combinations(["1", "2", "3", "4"], 3)) ==
  [["1", "2", "3"], ["1", "2", "4"], ["1", "3", "4"], ["2", "3", "4"]]
  @test collect(combinations(["1", "2", "3", "4"], 4)) == [["1", "2", "3", "4"]]
  @test collect(combinations(["1", "2", "3", "4"], 5)) == String[]
  @test collect(combinations(["1", "2", "3", "4"], 6)) == String[]

  v = collect(combinations(5, 3))
  @test issorted(v)
  @test allunique(v)

  @test collect(combinations(0, 0)) == [Int[]]
end
