@testset "ordered multiindices" begin
  I = Oscar.ordered_multi_index([1, 2, 5], 5)
  I1 = Oscar.ordered_multi_index([1], 5)
  I2 = Oscar.ordered_multi_index([2], 5)
  I5 = Oscar.ordered_multi_index([5], 5)
  J = Oscar.ordered_multi_index([3, 4], 5)
  K = Oscar.ordered_multi_index([2, 4], 5)
  @test "$(Oscar._wedge(I, J)[2])" == "0 < 1 < 2 < 3 < 4 < 5 <= 5"
  sign, ind = Oscar._wedge([I2, I5, I1])
  @test ind == I
  @test sign == 1
  sign, ind = Oscar._wedge([I2, I1, I5])
  @test ind == I
  @test sign == -1

  sign, ind = Oscar._wedge(J, K)
  @test sign == 0
  sign, ind = Oscar._wedge(I, K)
  @test sign == 0

  c = [ind for ind in Oscar.OrderedMultiIndexSet(3, 5)]
  @test length(c) == binomial(5, 3)

  @test all(x->c[Oscar.linear_index(x)] == x, Oscar.OrderedMultiIndexSet(3, 5))

  @test [Oscar.ordered_multi_index(k, 3, 5) for k in 1:binomial(5, 3)] == collect(Oscar.OrderedMultiIndexSet(3, 5))
end

@testset "multiindices of given degree" begin
  l = collect(Oscar.MultiIndicesOfDegree(3, 5)) # All non-negative multiindices I = [i1, i2, i3] 
                                              # with ∑ iₖ = 5
  @test length(l) == 21
  @test all(x->sum(x)==5, l)
  b = unique!(l)
  @test length(b) == length(l)

  l = collect(Oscar.MultiIndicesOfDegree(0, 5))
  @test l == Vector{Vector{Int}}()
  l = collect(Oscar.MultiIndicesOfDegree(7, 0))
  @test l == [[0 for i in 1:7]]
  l = collect(Oscar.MultiIndicesOfDegree(7, 1))
  @test length(l) == length(unique!(l)) == 7

  l = collect(Oscar.MultiIndicesOfDegree(1, 7))
  @test l == [[7]]
end

