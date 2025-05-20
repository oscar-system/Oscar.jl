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

  @testset "ranking/unranking" begin
    @test Oscar.combination(0,0,1) == Combination([])
    @test Oscar.combination(3,3,1) == Combination(collect(1:3))
    @test all(x->Oscar.combination(3,1,x) == Combination([x]), 1:3)

    n = 5
    k = 3

    C = combinations(n,k)
    c = collect(C)
    @test length(c) == binomial(n,k)
    @test all(x->c[Oscar.linear_index(x, n)] == x, combinations(n,k))
    @test all(i->C[i] == c[i], 1:length(c))

  end





  @testset "wedge products" begin
    n = 5
    k = 3

    I = Combination([1, 2, 5])
    I1 = Combination([1])
    I2 = Combination([2])
    I5 = Combination([5])
    J = Combination([3, 4])
    K = Combination([2, 4])

    sign, ind = Oscar._wedge(I, K)
    @test sign == 0
    sign, ind = Oscar._wedge(I, J)
    @test sign == 1
    ind == Combination(collect(1:5))
    sign, ind = Oscar._wedge(J, K)
    @test sign == 0

    sign, ind = Oscar._wedge([I2, I5, I1])
    @test ind == I
    @test sign == 1
    sign, ind = Oscar._wedge([I2, I1, I5])
    @test ind == I
    @test sign == -1
  end

end
