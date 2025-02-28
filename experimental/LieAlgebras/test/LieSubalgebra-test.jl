@testset "LieAlgebras.LieSubalgebra" begin
  @testset "Constructor and basic properties" begin
    let
      L = general_linear_lie_algebra(QQ, 3)
      b = sub(
        L, [basis(L, 1), basis(L, 2), basis(L, 3), basis(L, 5), basis(L, 6), basis(L, 9)]
      )

      n = sub(L, [basis(b, 2), basis(b, 3), basis(b, 5)])

      for S in [b, n]
        @test L == base_lie_algebra(S)
        @test length(gens(S)) == ngens(S)
        @test length(basis(S)) == dim(S)
        @test all(in(S), gens(S))
        @test all(in(S), basis(S))

        LS, emb = lie_algebra(S)
        @test dim(LS) == dim(S)
        @test domain(emb) == LS
        @test codomain(emb) == L
        @test image(emb) == S
      end

      @test dim(b) == 6
      @test dim(n) == 3

      @test n == sub(L, [basis(b, 2), basis(b, 3), basis(b, 5)])
      @test n == sub(L, [basis(b, 2) + basis(b, 3), basis(b, 3), basis(b, 5)])
      @test n ==
        sub(L, [basis(b, 2) + basis(b, 3), basis(b, 3), basis(b, 5)]; is_basis=true)
    end
  end

  @testset "Hum72, Exercise 2.3" begin
    @testset for n in 2:4, F in [QQ, GF(2)]
      L = general_linear_lie_algebra(F, n)

      b = sub(L, [basis(L, (i - 1) * n + j) for i in 1:n, j in 1:n if i <= j])
      d = sub(L, [basis(L, (i - 1) * n + j) for i in 1:n, j in 1:n if i == j])
      n = sub(L, [basis(L, (i - 1) * n + j) for i in 1:n, j in 1:n if i < j])

      @test is_self_normalizing(b)
      @test is_self_normalizing(d)
      @test !is_self_normalizing(n)
      @test normalizer(n) == b
    end
  end

  @testset "#4676" begin
    L = special_linear_lie_algebra(QQ, 2)
    e,f,h = basis(L)
    sub1 = sub(L, e)
    sub2 = sub(L, f)
    sub3 = bracket(sub1, sub2)
    @test sub3 isa LieSubalgebra
    @test dim(sub3) < dim(L)
  end
end
